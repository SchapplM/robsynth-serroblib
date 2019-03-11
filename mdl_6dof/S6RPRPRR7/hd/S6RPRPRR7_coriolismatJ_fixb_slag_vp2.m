% Calculate matrix of centrifugal and coriolis load on the joints for
% S6RPRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 03:56
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Cq = S6RPRPRR7_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR7_coriolismatJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR7_coriolismatJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR7_coriolismatJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR7_coriolismatJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR7_coriolismatJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR7_coriolismatJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:55:03
% EndTime: 2019-03-09 03:55:17
% DurationCPUTime: 8.18s
% Computational Cost: add. (22308->506), mult. (39460->671), div. (0->0), fcn. (46193->8), ass. (0->283)
t283 = sin(pkin(10));
t286 = sin(qJ(3));
t288 = cos(qJ(3));
t413 = cos(pkin(10));
t257 = t283 * t286 - t288 * t413;
t258 = -t283 * t288 - t286 * t413;
t285 = sin(qJ(5));
t455 = cos(qJ(5));
t492 = -t455 * t257 + t258 * t285;
t514 = -t492 / 0.2e1;
t540 = 0.2e1 * t514;
t284 = sin(qJ(6));
t458 = t284 / 0.2e1;
t287 = cos(qJ(6));
t456 = t287 / 0.2e1;
t225 = t285 * t257 + t455 * t258;
t452 = pkin(5) * t225;
t281 = t284 ^ 2;
t282 = t287 ^ 2;
t370 = t281 + t282;
t503 = t370 * t492;
t539 = pkin(9) * t503 + t452;
t274 = t286 * pkin(3) + qJ(2);
t238 = -pkin(4) * t258 + t274;
t459 = -t284 / 0.2e1;
t468 = -t225 / 0.2e1;
t277 = Ifges(7,5) * t287;
t442 = Ifges(7,6) * t284;
t490 = t277 - t442;
t278 = Ifges(7,4) * t287;
t489 = -Ifges(7,2) * t284 + t278;
t498 = Ifges(7,6) * t492 + t225 * t489;
t445 = Ifges(7,4) * t284;
t268 = Ifges(7,1) * t287 - t445;
t499 = Ifges(7,5) * t492 + t225 * t268;
t513 = t492 / 0.2e1;
t538 = Ifges(6,4) * t540 + t456 * t499 + t459 * t498 + t238 * mrSges(6,1) + t490 * t513 + (Ifges(6,2) + Ifges(7,3)) * t468;
t509 = t492 * mrSges(6,2);
t529 = t225 * mrSges(6,1);
t537 = t529 - t509;
t422 = t287 * mrSges(7,2);
t430 = t284 * mrSges(7,1);
t262 = t422 + t430;
t504 = t262 * t492;
t179 = -t504 / 0.2e1;
t318 = t422 / 0.2e1 + t430 / 0.2e1;
t313 = t318 * t492;
t523 = t179 - t313;
t536 = t523 * qJD(6);
t535 = t498 * t456 + t458 * t499;
t423 = t287 * mrSges(7,1);
t429 = t284 * mrSges(7,2);
t334 = t423 - t429;
t527 = t225 * t334;
t534 = t527 / 0.2e1 + t529 / 0.2e1;
t443 = Ifges(6,6) * t492;
t265 = Ifges(7,2) * t287 + t445;
t382 = t284 * t265;
t267 = Ifges(7,1) * t284 + t278;
t532 = t267 * t456;
t491 = -t382 / 0.2e1 + t532;
t533 = (Ifges(6,5) + t491) * t225 - t443 + t535;
t531 = t504 / 0.2e1;
t475 = pkin(1) + pkin(7);
t369 = qJ(4) + t475;
t340 = t283 * t369;
t251 = t288 * t340;
t322 = t369 * t413;
t230 = -t286 * t322 - t251;
t449 = t258 * pkin(8);
t184 = t230 + t449;
t316 = t288 * t322;
t228 = t286 * t340 - t316;
t450 = t257 * pkin(8);
t299 = t228 + t450;
t116 = t184 * t285 - t299 * t455;
t505 = t262 * t225;
t528 = t116 * t505;
t524 = t492 * t225;
t453 = pkin(5) * t492;
t156 = -pkin(9) * t225 + t453;
t522 = t514 + t513;
t395 = t492 * t284;
t148 = mrSges(7,2) * t225 - mrSges(7,3) * t395;
t377 = t287 * t148;
t394 = t492 * t287;
t151 = -mrSges(7,1) * t225 - mrSges(7,3) * t394;
t383 = t284 * t151;
t317 = t377 / 0.2e1 - t383 / 0.2e1;
t118 = t184 * t455 + t285 * t299;
t123 = -pkin(9) * t492 + t238 - t452;
t53 = -t118 * t284 + t123 * t287;
t54 = t118 * t287 + t123 * t284;
t329 = t284 * t53 - t287 * t54;
t478 = m(7) / 0.2e1;
t520 = -t329 * t478 + t317;
t339 = t286 * t369;
t231 = t283 * t339 - t316;
t185 = t231 + t450;
t229 = t339 * t413 + t251;
t307 = t229 - t449;
t119 = t185 * t455 + t285 * t307;
t117 = t185 * t285 - t307 * t455;
t324 = -t116 * t225 - t117 * t492;
t519 = -t118 * t492 + t119 * t225 - t324;
t517 = t225 / 0.2e1;
t516 = t225 / 0.4e1;
t515 = -t262 / 0.2e1;
t511 = mrSges(7,1) * t492;
t510 = mrSges(7,2) * t492;
t434 = t225 * mrSges(6,3);
t436 = t492 * mrSges(6,3);
t493 = (t277 / 0.2e1 - t442 / 0.2e1) * t225;
t502 = -t257 * t283 + t413 * t258;
t280 = t288 * pkin(3);
t239 = -pkin(4) * t257 + t280;
t124 = t239 + t156;
t56 = t119 * t287 + t124 * t284;
t417 = t56 * t287;
t55 = -t119 * t284 + t124 * t287;
t418 = t55 * t284;
t328 = t417 - t418;
t501 = t286 * mrSges(4,1) + t288 * mrSges(4,2);
t500 = -t268 / 0.4e1 + t265 / 0.4e1;
t497 = -Ifges(7,6) / 0.2e1;
t273 = pkin(3) * t413 + pkin(4);
t454 = pkin(3) * t283;
t243 = t285 * t273 + t454 * t455;
t241 = pkin(9) + t243;
t342 = t370 * t241;
t376 = t287 * t151;
t386 = t284 * t148;
t488 = -t376 / 0.2e1 - t386 / 0.2e1;
t428 = t284 * mrSges(7,3);
t147 = -t225 * t428 - t510;
t378 = t287 * t147;
t433 = t225 * mrSges(7,3);
t150 = -t287 * t433 + t511;
t384 = t284 * t150;
t487 = t378 / 0.2e1 - t384 / 0.2e1;
t398 = t225 * t284;
t146 = -mrSges(7,3) * t398 - t510;
t379 = t287 * t146;
t375 = t287 * t225;
t149 = -mrSges(7,3) * t375 + t511;
t385 = t284 * t149;
t486 = t379 / 0.2e1 - t385 / 0.2e1;
t101 = -Ifges(7,5) * t225 + t268 * t492;
t380 = t287 * t101;
t98 = -Ifges(7,6) * t225 + t489 * t492;
t426 = t284 * t98;
t485 = Ifges(6,1) * t513 + Ifges(6,4) * t225 + t380 / 0.2e1 - t426 / 0.2e1;
t484 = t116 * t515 + t490 * t516;
t448 = mrSges(7,3) * t492;
t483 = t370 * t448;
t482 = 0.2e1 * m(7);
t481 = -m(5) / 0.2e1;
t480 = m(6) / 0.2e1;
t479 = -m(7) / 0.2e1;
t477 = -pkin(5) / 0.2e1;
t476 = m(5) * pkin(3);
t371 = t229 + t230;
t372 = t228 - t231;
t6 = t504 * t468 - (-t436 / 0.2e1 + t487) * t225 + (-t225 * t328 + t324) * t478 - t519 * t480 + m(5) * (-t257 * t371 + t258 * t372) / 0.2e1 + (-t434 / 0.2e1 - t505 / 0.2e1 + t520) * t492;
t63 = -t116 * t287 + t156 * t284;
t415 = t63 * t287;
t62 = t116 * t284 + t156 * t287;
t416 = t62 * t284;
t327 = t415 - t416;
t9 = -((t116 + t327) * t478 + t531 + t486) * t225 - ((t118 + t329) * t478 + t505 / 0.2e1 - t317) * t492;
t474 = t6 * qJD(3) + t9 * qJD(5);
t467 = -t225 / 0.4e1;
t242 = t273 * t455 - t285 * t454;
t240 = -pkin(5) - t242;
t463 = -t240 / 0.2e1;
t462 = -t241 / 0.2e1;
t461 = t241 / 0.2e1;
t460 = -t242 / 0.2e1;
t457 = -t287 / 0.2e1;
t451 = pkin(5) * t262;
t441 = t116 * mrSges(6,2);
t440 = t117 * mrSges(6,1);
t439 = t118 * mrSges(6,1);
t438 = t119 * mrSges(6,2);
t437 = t492 * mrSges(6,1);
t427 = t284 * t54;
t140 = t334 * t492;
t144 = t265 * t492;
t145 = t492 * t267;
t263 = Ifges(7,5) * t284 + Ifges(7,6) * t287;
t7 = t116 * t140 + t53 * t148 - t54 * t151 + (t329 * mrSges(7,3) - t145 * t456 + t263 * t517 + t98 * t457 + (t101 - t144) * t459) * t492;
t414 = t7 * qJD(1);
t348 = t257 ^ 2 + t258 ^ 2;
t409 = t116 * t492;
t14 = (t504 + t436) * t492 + t348 * mrSges(5,3) + (t377 - t383 + t434) * t225 + m(7) * (-t225 * t329 + t409) + m(6) * (t118 * t225 + t409) + m(5) * (t228 * t257 + t230 * t258);
t412 = qJD(1) * t14;
t346 = -t258 * mrSges(5,1) - t257 * mrSges(5,2);
t304 = -t346 - t501 + t537;
t26 = t386 + t376 + mrSges(3,3) + (m(4) + m(3)) * qJ(2) + m(7) * (t287 * t53 + t427) + m(6) * t238 + m(5) * t274 - t304;
t411 = qJD(1) * t26;
t410 = t116 * t117;
t407 = t117 * t334;
t406 = t118 * t334;
t220 = t225 * mrSges(6,2);
t247 = t258 * mrSges(5,2);
t360 = t334 * t513;
t368 = t476 / 0.2e1;
t291 = (t225 * t243 - t242 * t492) * t480 + (t225 * t342 + t240 * t492) * t478 - t360 + (t257 * t413 + t258 * t283) * t368 + t370 * t433 / 0.2e1;
t295 = t239 * t480 + (t284 * t56 + t287 * t55) * t478 + t147 * t458 + t150 * t456 + t288 * t368;
t17 = t257 * mrSges(5,1) - t220 - t247 + t291 - t295 - t437;
t405 = t17 * qJD(1);
t349 = -t282 / 0.2e1 - t281 / 0.2e1;
t303 = t349 * t448 + t488;
t296 = t140 * t514 - t225 * t303;
t319 = -t429 / 0.2e1 + t423 / 0.2e1;
t22 = t296 - t319;
t404 = t22 * qJD(1);
t403 = t492 ^ 2;
t396 = t492 * t263;
t336 = mrSges(7,3) * t349;
t344 = t370 * t225;
t294 = -t225 * t336 - t220 / 0.2e1 + (pkin(9) * t344 - t453) * t478 - t360;
t297 = (t284 * t63 + t287 * t62) * t479 + mrSges(6,2) * t468 + t146 * t459 + t149 * t457;
t24 = mrSges(6,1) * t540 + t294 + t297;
t393 = t24 * qJD(1);
t392 = t240 * t140;
t391 = t240 * t262;
t314 = t318 * t225;
t27 = -t314 - t317;
t388 = t27 * qJD(1);
t293 = (-t225 * t344 - t403) * t478 + (-t225 ^ 2 - t403) * t480 + t348 * t481;
t315 = t481 - m(6) / 0.2e1 + t370 * t479;
t36 = t293 + t315;
t373 = t36 * qJD(1);
t367 = mrSges(7,3) * t415;
t366 = t505 * t477;
t365 = Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t362 = -t428 / 0.2e1;
t361 = t334 / 0.2e1 + mrSges(6,1) / 0.2e1;
t350 = t516 + t467;
t341 = t370 * t242;
t337 = t268 * t458 + t456 * t489 + t491;
t1 = t528 + t117 * t504 + t54 * t147 + t56 * t148 + t53 * t150 + t55 * t151 - t239 * t529 + t238 * t220 + t274 * t247 + (t468 * t490 + t485) * t225 + (-t274 * mrSges(5,1) + t371 * mrSges(5,3) - Ifges(5,4) * t257) * t257 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t288 + pkin(3) * t346) * t288 + m(5) * (t228 * t229 + t230 * t231 + t274 * t280) + m(6) * (t118 * t119 + t238 * t239 + t410) + m(7) * (t53 * t55 + t54 * t56 + t410) + t519 * mrSges(6,3) + (-mrSges(4,2) * qJ(2) + Ifges(4,4) * t286 + (-Ifges(4,1) + Ifges(4,2)) * t288) * t286 + (Ifges(5,4) * t258 + (-Ifges(5,1) + Ifges(5,2)) * t257 - t372 * mrSges(5,3)) * t258 + (t239 * mrSges(6,2) + Ifges(6,1) * t517 - t225 * t365 + t538) * t492;
t331 = t1 * qJD(1) + t6 * qJD(2);
t4 = m(7) * (t116 * t118 + t53 * t62 + t54 * t63) + t63 * t148 + t54 * t146 + t62 * t151 + t53 * t149 + t118 * t504 + t528 + t538 * t492 - (-t238 * mrSges(6,2) + t493 + (-Ifges(6,1) / 0.2e1 + t365) * t492 - t485) * t225;
t330 = t4 * qJD(1) + t9 * qJD(2);
t37 = m(7) * (-t225 * t503 + t524);
t326 = t6 * qJD(1) + t37 * qJD(2);
t42 = m(7) * (0.1e1 - t370) * t524;
t325 = t9 * qJD(1) + t42 * qJD(2);
t323 = t116 * t243 + t118 * t240;
t298 = m(7) * (-(t240 + t341) * t225 - (t243 - t342) * t492);
t300 = t478 * t539 + t534;
t19 = -t361 * t225 + t522 * mrSges(6,2) - t298 / 0.2e1 + t300;
t302 = t439 / 0.2e1 + t406 / 0.2e1 + t505 * t463 + t243 * t179;
t2 = t366 + t323 * t479 + t522 * Ifges(6,6) + (t517 + t468) * Ifges(6,5) + (-t119 / 0.2e1 - t116 / 0.2e1) * mrSges(6,2) + (m(7) * t477 - t361) * t117 + (t148 * t460 + t146 * t462 + pkin(9) * t147 / 0.2e1 + t350 * t267 + (-t63 / 0.2e1 + t56 / 0.2e1) * mrSges(7,3) + (pkin(9) * t56 / 0.4e1 - t241 * t63 / 0.4e1 - t242 * t54 / 0.4e1) * t482) * t287 + (-pkin(9) * t150 / 0.2e1 + t242 * t151 / 0.2e1 + t149 * t461 - t350 * t265 + (t62 / 0.2e1 - t55 / 0.2e1) * mrSges(7,3) + (-pkin(9) * t55 / 0.4e1 + t241 * t62 / 0.4e1 + t242 * t53 / 0.4e1) * t482) * t284 + t302;
t301 = -t242 * mrSges(6,2) + mrSges(7,3) * t341 + (-t334 - mrSges(6,1)) * t243;
t44 = m(7) * (t240 * t243 + t241 * t341) + t301;
t312 = -t2 * qJD(1) - t19 * qJD(2) + t44 * qJD(3);
t309 = Ifges(7,3) * t513 + t55 * mrSges(7,1) / 0.2e1 - t56 * mrSges(7,2) / 0.2e1;
t10 = -t392 / 0.2e1 + (Ifges(7,5) * t517 + t144 / 0.4e1 - t101 / 0.4e1 + t151 * t461) * t287 + (t225 * t497 + t98 / 0.4e1 + t145 / 0.4e1 + t148 * t461) * t284 + (t500 * t287 + (t267 / 0.4e1 + t489 / 0.4e1) * t284 - t241 * t336) * t492 + t309 + t484;
t125 = t337 + t391;
t68 = t531 - t313;
t310 = -t10 * qJD(1) - t68 * qJD(2) + t125 * qJD(3);
t308 = Ifges(7,3) * t514 - t62 * mrSges(7,1) / 0.2e1 + t63 * mrSges(7,2) / 0.2e1;
t305 = t54 * t362 - t284 * t145 / 0.4e1 - t287 * t144 / 0.4e1 + mrSges(7,3) * t427 / 0.2e1 - t426 / 0.4e1 + t380 / 0.4e1 - t484 - t500 * t394 - (t489 + t267) * t395 / 0.4e1;
t292 = pkin(9) * t303 + t140 * t477 + t305;
t13 = t292 + t308 - t493;
t161 = t337 - t451;
t66 = (t463 + pkin(5) / 0.2e1) * t262 + (mrSges(7,2) * t460 - t267 / 0.2e1 - t489 / 0.2e1) * t287 + (mrSges(7,1) * t460 - t268 / 0.2e1 + t265 / 0.2e1) * t284;
t69 = (t515 + t318) * t492;
t306 = t13 * qJD(1) + t69 * qJD(2) - t66 * qJD(3) + t161 * qJD(5);
t67 = t391 / 0.2e1 - t451 / 0.2e1 - t318 * t242 + t337;
t35 = t293 - t315;
t28 = -t314 + t317;
t25 = mrSges(6,1) * t513 - t437 / 0.2e1 + t294 - t297;
t23 = t296 + t319;
t21 = t291 + t295;
t20 = t298 / 0.2e1 - t509 + t300 + (-t349 * t492 + t370 * t513) * mrSges(7,3) + t534;
t12 = Ifges(7,5) * t375 / 0.2e1 + t398 * t497 + t292 - t308;
t11 = t493 + t305 + t392 / 0.2e1 + t309 + t462 * t483 + t488 * t241;
t3 = -t438 / 0.2e1 - t302 - t443 / 0.2e1 + t367 / 0.2e1 + t396 / 0.4e1 + t62 * t362 + (t417 / 0.2e1 - t418 / 0.2e1) * mrSges(7,3) - t407 / 0.2e1 - t440 / 0.2e1 + t441 / 0.2e1 + t366 + (t263 / 0.4e1 - Ifges(6,6) / 0.2e1) * t492 + t382 * t467 + Ifges(6,5) * t517 + (-pkin(5) * t117 + t323) * t478 + t520 * t242 + (t327 * t478 + t486) * t241 + (t328 * t478 + t487) * pkin(9) + (Ifges(6,5) / 0.2e1 + t532 - t382 / 0.4e1) * t225 + t535;
t5 = [qJD(2) * t26 + qJD(3) * t1 + qJD(4) * t14 + qJD(5) * t4 + qJD(6) * t7, qJD(4) * t35 + qJD(6) * t23 + t411 + t474, t21 * qJD(4) + t3 * qJD(5) + t11 * qJD(6) + t331 + (-t438 - t502 * mrSges(5,3) * pkin(3) + m(6) * (-t117 * t242 + t119 * t243) + t501 * t475 + t328 * mrSges(7,3) - Ifges(4,6) * t288 - Ifges(4,5) * t286 + Ifges(5,5) * t258 + Ifges(5,6) * t257 + t229 * mrSges(5,1) - t231 * mrSges(5,2) - t243 * t436 - t242 * t434 - t407 + t533 - t440 + t263 * t513 + (m(7) * t117 + t505) * t240 + (m(7) * t328 + t378 - t384) * t241 + (t229 * t413 + t231 * t283) * t476) * qJD(3), qJD(2) * t35 + qJD(3) * t21 + qJD(5) * t25 + qJD(6) * t28 + t412, t3 * qJD(3) + t25 * qJD(4) + t12 * qJD(6) + t330 + (-t406 + t367 - mrSges(7,3) * t416 + t396 / 0.2e1 + t441 - t439 + (-m(7) * t118 - t505) * pkin(5) + (m(7) * t327 + t379 - t385) * pkin(9) + t533) * qJD(5), t414 + t23 * qJD(2) + t11 * qJD(3) + t28 * qJD(4) + t12 * qJD(5) + (-mrSges(7,1) * t54 - mrSges(7,2) * t53 - t396) * qJD(6); qJD(4) * t36 + qJD(6) * t22 - t411 + t474, qJD(3) * t37 + qJD(5) * t42 (m(7) * (-t225 * t240 + t342 * t492) + t527 + m(6) * (t225 * t242 + t243 * t492) + t502 * t476 + t304 + t483) * qJD(3) + t20 * qJD(5) + t536 + t326, t373, t20 * qJD(3) + (m(7) * t539 + mrSges(7,3) * t503 + t527 + t537) * qJD(5) + t536 + t325, qJD(6) * t527 + t404 + (qJD(3) + qJD(5)) * t523; qJD(4) * t17 - qJD(5) * t2 - qJD(6) * t10 - t331, -qJD(5) * t19 - qJD(6) * t68 - t326, qJD(5) * t44 + qJD(6) * t125, t405 (m(7) * (-pkin(5) * t243 + pkin(9) * t341) + t301) * qJD(5) + t67 * qJD(6) + t312, t67 * qJD(5) + (-t241 * t334 + t490) * qJD(6) + t310; -qJD(2) * t36 - qJD(3) * t17 - qJD(5) * t24 - qJD(6) * t27 - t412, -t373, -t405, 0, -t393, -qJD(6) * t262 - t388; qJD(3) * t2 + qJD(4) * t24 + qJD(6) * t13 - t330, qJD(3) * t19 + qJD(6) * t69 - t325, -qJD(6) * t66 - t312, t393, t161 * qJD(6) (-pkin(9) * t334 + t490) * qJD(6) + t306; -qJD(2) * t22 + qJD(3) * t10 + qJD(4) * t27 - qJD(5) * t13 - t414, qJD(3) * t68 - qJD(5) * t69 - t404, qJD(5) * t66 - t310, t388, -t306, 0;];
Cq  = t5;
