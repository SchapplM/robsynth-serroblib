% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:39:11
% EndTime: 2019-03-09 05:40:10
% DurationCPUTime: 28.97s
% Computational Cost: add. (30581->826), mult. (101482->1176), div. (0->0), fcn. (87253->14), ass. (0->368)
t312 = sin(pkin(12));
t314 = sin(pkin(6));
t316 = cos(pkin(12));
t320 = sin(qJ(3));
t323 = cos(qJ(3));
t313 = sin(pkin(7));
t403 = cos(pkin(6));
t356 = t403 * t313;
t317 = cos(pkin(7));
t387 = t317 * t323;
t324 = t323 * t356 + t314 * (-t312 * t320 + t316 * t387);
t246 = t324 * qJD(1);
t244 = qJD(4) - t246;
t388 = t317 * t320;
t330 = (-t312 * t388 + t316 * t323) * t314;
t269 = qJD(1) * t330;
t376 = qJD(3) * t323;
t521 = t313 * t376 - t269;
t364 = pkin(1) * t403;
t352 = qJD(1) * t364;
t303 = t316 * t352;
t396 = t312 * t314;
t326 = t403 * pkin(2) + (-pkin(9) * t317 - qJ(2)) * t396;
t245 = qJD(1) * t326 + t303;
t271 = (-pkin(9) * t312 * t313 - pkin(2) * t316 - pkin(1)) * t314;
t264 = qJD(1) * t271 + qJD(2);
t340 = t245 * t317 + t264 * t313;
t392 = t314 * t316;
t305 = qJ(2) * t392;
t278 = qJD(1) * t305 + t312 * t352;
t329 = (t317 * t392 + t356) * pkin(9);
t238 = qJD(1) * t329 + t278;
t386 = t323 * t238;
t176 = t320 * t340 + t386;
t319 = sin(qJ(4));
t322 = cos(qJ(4));
t520 = -qJD(5) * t319 - t176 + t244 * (pkin(4) * t319 - qJ(5) * t322);
t394 = t313 * t320;
t287 = -t322 * t317 + t319 * t394;
t379 = qJD(1) * t314;
t363 = t312 * t379;
t354 = t313 * t363;
t474 = -qJD(4) * t287 - t319 * t354 + t521 * t322;
t331 = (t312 * t387 + t316 * t320) * t314;
t268 = qJD(1) * t331;
t377 = qJD(3) * t320;
t519 = t313 * t377 - t268;
t247 = t324 * qJD(3);
t236 = qJD(1) * t247;
t254 = t320 * t356 + (t312 * t323 + t316 * t388) * t314;
t249 = t254 * qJD(1);
t282 = -t313 * t392 + t317 * t403;
t325 = qJD(1) * t282 + qJD(3);
t215 = t319 * t249 - t322 * t325;
t375 = qJD(4) * t215;
t180 = t322 * t236 - t375;
t446 = t180 / 0.2e1;
t518 = Ifges(5,4) * t446;
t227 = t320 * t238;
t175 = t340 * t323 - t227;
t207 = pkin(3) * t249 - pkin(10) * t246;
t119 = t322 * t175 + t319 * t207;
t106 = qJ(5) * t249 + t119;
t311 = sin(pkin(13));
t315 = cos(pkin(13));
t374 = qJD(4) * t319;
t372 = pkin(10) * t374;
t486 = t520 * t315 + (t106 + t372) * t311;
t517 = -t315 * t106 + t311 * t520;
t216 = t322 * t249 + t319 * t325;
t181 = qJD(4) * t216 + t319 * t236;
t445 = -t181 / 0.2e1;
t248 = t254 * qJD(3);
t237 = qJD(1) * t248;
t435 = t237 / 0.2e1;
t390 = t315 * t322;
t202 = t246 * t390 + t249 * t311;
t400 = t246 * t319;
t516 = -pkin(5) * t400 + pkin(11) * t202 + (pkin(5) * t319 - pkin(11) * t390) * qJD(4) + t486;
t397 = t311 * t322;
t201 = -t246 * t397 + t249 * t315;
t391 = t315 * t319;
t515 = -pkin(11) * t201 + (-pkin(10) * t391 - pkin(11) * t397) * qJD(4) + t517;
t476 = -t474 * t311 + t519 * t315;
t475 = t519 * t311 + t474 * t315;
t208 = -t245 * t313 + t317 * t264;
t148 = -pkin(3) * t246 - pkin(10) * t249 + t208;
t151 = pkin(10) * t325 + t245 * t388 + t264 * t394 + t386;
t85 = t148 * t319 + t151 * t322;
t80 = qJ(5) * t244 + t85;
t393 = t313 * t323;
t150 = -pkin(3) * t325 - t245 * t387 - t264 * t393 + t227;
t96 = t215 * pkin(4) - t216 * qJ(5) + t150;
t35 = -t311 * t80 + t315 * t96;
t36 = t311 * t96 + t315 * t80;
t185 = t216 * t315 + t244 * t311;
t25 = pkin(5) * t215 - pkin(11) * t185 + t35;
t355 = -t216 * t311 + t315 * t244;
t27 = pkin(11) * t355 + t36;
t318 = sin(qJ(6));
t321 = cos(qJ(6));
t5 = t25 * t321 - t27 * t318;
t6 = t25 * t318 + t27 * t321;
t468 = -t5 * mrSges(7,1) + t6 * mrSges(7,2);
t514 = -t35 * mrSges(6,1) + t36 * mrSges(6,2) + t468;
t327 = qJD(2) * t330;
t143 = qJD(1) * t327 + qJD(3) * t175;
t378 = qJD(2) * t314;
t362 = t312 * t378;
t351 = qJD(1) * t362;
t336 = t313 * t351;
t195 = pkin(3) * t237 - pkin(10) * t236 + t336;
t373 = qJD(4) * t322;
t42 = t322 * t143 + t148 * t373 - t151 * t374 + t319 * t195;
t33 = qJ(5) * t237 + qJD(5) * t244 + t42;
t328 = qJD(2) * t331;
t144 = qJD(1) * t328 + qJD(3) * t176;
t67 = t181 * pkin(4) - t180 * qJ(5) - t216 * qJD(5) + t144;
t21 = -t311 * t33 + t315 * t67;
t22 = t311 * t67 + t315 * t33;
t444 = t181 / 0.2e1;
t138 = t180 * t315 + t237 * t311;
t447 = t138 / 0.2e1;
t137 = -t180 * t311 + t237 * t315;
t448 = t137 / 0.2e1;
t113 = t185 * t321 + t318 * t355;
t48 = -qJD(6) * t113 + t137 * t321 - t138 * t318;
t461 = t48 / 0.2e1;
t507 = -t185 * t318 + t321 * t355;
t47 = qJD(6) * t507 + t137 * t318 + t138 * t321;
t462 = t47 / 0.2e1;
t10 = pkin(11) * t137 + t22;
t7 = pkin(5) * t181 - pkin(11) * t138 + t21;
t1 = qJD(6) * t5 + t10 * t321 + t318 * t7;
t2 = -qJD(6) * t6 - t10 * t318 + t321 * t7;
t469 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t513 = t469 - 0.2e1 * Ifges(5,2) * t445 - 0.2e1 * Ifges(5,6) * t435 - t22 * mrSges(6,2) - (-Ifges(7,3) - Ifges(6,3)) * t444 - t518 + Ifges(6,5) * t447 + Ifges(7,5) * t462 + Ifges(6,6) * t448 + Ifges(7,6) * t461 + t21 * mrSges(6,1);
t409 = t244 * Ifges(5,6);
t423 = Ifges(5,4) * t216;
t427 = t85 * mrSges(5,3);
t512 = t427 + t409 / 0.2e1 + t423 / 0.2e1 - t150 * mrSges(5,1);
t84 = t148 * t322 - t319 * t151;
t511 = t84 * mrSges(5,3);
t509 = t150 * mrSges(5,2);
t380 = t312 * t364 + t305;
t250 = t329 + t380;
t307 = t316 * t364;
t255 = t307 + t326;
t339 = t255 * t317 + t271 * t313;
t188 = -t320 * t250 + t339 * t323;
t126 = -t215 * Ifges(5,2) + t409 + t423;
t412 = t185 * Ifges(6,5);
t413 = t355 * Ifges(6,6);
t508 = t126 / 0.2e1 - t413 / 0.2e1 - t412 / 0.2e1 + t514;
t505 = Ifges(5,1) * t446 + Ifges(5,5) * t435;
t456 = Ifges(5,4) * t445 + t505;
t214 = qJD(6) + t215;
t439 = t214 / 0.2e1;
t452 = t113 / 0.2e1;
t454 = t507 / 0.2e1;
t411 = t214 * Ifges(7,3);
t414 = t113 * Ifges(7,5);
t415 = t507 * Ifges(7,6);
t53 = t411 + t414 + t415;
t91 = t215 * Ifges(6,3) + t412 + t413;
t495 = t91 + t53;
t504 = -t126 / 0.2e1 + Ifges(7,5) * t452 + Ifges(7,6) * t454 + Ifges(7,3) * t439 + t495 / 0.2e1;
t502 = t375 / 0.2e1;
t499 = t84 * mrSges(5,1);
t498 = t85 * mrSges(5,2);
t11 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t181;
t496 = t138 * Ifges(6,5) + t137 * Ifges(6,6) + t181 * Ifges(6,3) + t11;
t494 = Ifges(4,5) * t236;
t493 = Ifges(4,6) * t237;
t492 = t143 * mrSges(4,2);
t141 = mrSges(5,1) * t237 - mrSges(5,3) * t180;
t82 = -t137 * mrSges(6,1) + t138 * mrSges(6,2);
t491 = -t141 + t82;
t297 = -pkin(4) * t322 - qJ(5) * t319 - pkin(3);
t291 = t315 * t297;
t251 = -pkin(11) * t391 + t291 + (-pkin(10) * t311 - pkin(5)) * t322;
t273 = pkin(10) * t390 + t311 * t297;
t398 = t311 * t319;
t263 = -pkin(11) * t398 + t273;
t210 = t251 * t318 + t263 * t321;
t490 = -qJD(6) * t210 - t515 * t318 + t516 * t321;
t209 = t251 * t321 - t263 * t318;
t489 = qJD(6) * t209 + t516 * t318 + t515 * t321;
t425 = pkin(11) + qJ(5);
t298 = t425 * t311;
t299 = t425 * t315;
t265 = -t298 * t321 - t299 * t318;
t337 = t311 * t318 - t315 * t321;
t164 = pkin(4) * t216 + qJ(5) * t215;
t68 = t315 * t164 - t311 * t84;
t39 = pkin(11) * t215 * t315 + pkin(5) * t216 + t68;
t401 = t215 * t311;
t69 = t311 * t164 + t315 * t84;
t56 = pkin(11) * t401 + t69;
t488 = -qJD(5) * t337 + qJD(6) * t265 - t318 * t39 - t321 * t56;
t266 = -t298 * t318 + t299 * t321;
t294 = t311 * t321 + t315 * t318;
t487 = -qJD(5) * t294 - qJD(6) * t266 + t318 * t56 - t321 * t39;
t485 = -t315 * t372 + t517;
t118 = -t319 * t175 + t207 * t322;
t107 = -pkin(4) * t249 - t118;
t365 = pkin(5) * t311 + pkin(10);
t484 = pkin(5) * t201 + t365 * t373 - t107;
t128 = t201 * t321 - t202 * t318;
t283 = t337 * qJD(6);
t231 = t283 * t319 - t294 * t373;
t480 = t128 - t231;
t129 = t201 * t318 + t202 * t321;
t284 = t294 * qJD(6);
t230 = -t284 * t319 - t337 * t373;
t479 = t129 - t230;
t288 = t317 * t319 + t322 * t394;
t258 = -t311 * t288 - t315 * t393;
t259 = t315 * t288 - t311 * t393;
t211 = t258 * t321 - t259 * t318;
t478 = qJD(6) * t211 + t476 * t318 + t475 * t321;
t212 = t258 * t318 + t259 * t321;
t477 = -qJD(6) * t212 - t475 * t318 + t476 * t321;
t407 = t249 * mrSges(4,3);
t383 = -mrSges(4,1) * t325 + mrSges(5,1) * t215 + mrSges(5,2) * t216 + t407;
t473 = qJD(4) * t288 + t521 * t319 + t322 * t354;
t43 = -t319 * t143 - t148 * t374 - t151 * t373 + t195 * t322;
t37 = -pkin(4) * t237 - t43;
t79 = -pkin(4) * t244 + qJD(5) - t84;
t472 = t319 * t37 + t79 * t373;
t471 = -t319 * t43 + t322 * t42;
t92 = t185 * Ifges(6,4) + Ifges(6,2) * t355 + Ifges(6,6) * t215;
t457 = -t92 / 0.2e1;
t470 = t311 * t457 - t511;
t467 = -m(5) * t150 - t383;
t466 = t43 * mrSges(5,1) - t42 * mrSges(5,2);
t464 = Ifges(7,4) * t462 + Ifges(7,2) * t461 + Ifges(7,6) * t444;
t463 = Ifges(7,1) * t462 + Ifges(7,4) * t461 + Ifges(7,5) * t444;
t418 = Ifges(7,4) * t113;
t54 = Ifges(7,2) * t507 + Ifges(7,6) * t214 + t418;
t460 = t54 / 0.2e1;
t109 = Ifges(7,4) * t507;
t55 = Ifges(7,1) * t113 + Ifges(7,5) * t214 + t109;
t459 = t55 / 0.2e1;
t458 = Ifges(6,1) * t447 + Ifges(6,4) * t448 + Ifges(6,5) * t444;
t455 = -t507 / 0.2e1;
t453 = -t113 / 0.2e1;
t213 = Ifges(5,4) * t215;
t410 = t244 * Ifges(5,5);
t127 = t216 * Ifges(5,1) - t213 + t410;
t449 = t127 / 0.2e1;
t443 = t355 / 0.2e1;
t442 = t185 / 0.2e1;
t440 = -t214 / 0.2e1;
t438 = -t215 / 0.2e1;
t437 = t215 / 0.2e1;
t432 = t249 / 0.2e1;
t431 = t315 / 0.2e1;
t160 = qJD(3) * t188 + t327;
t217 = -t255 * t313 + t317 * t271;
t167 = -pkin(3) * t324 - pkin(10) * t254 + t217;
t242 = t323 * t250;
t189 = t255 * t388 + t271 * t394 + t242;
t174 = pkin(10) * t282 + t189;
t353 = t313 * t362;
t199 = pkin(3) * t248 - pkin(10) * t247 + t353;
t59 = t322 * t160 + t167 * t373 - t174 * t374 + t319 * t199;
t41 = qJ(5) * t248 - qJD(5) * t324 + t59;
t161 = t328 + (t320 * t339 + t242) * qJD(3);
t221 = t254 * t322 + t282 * t319;
t193 = qJD(4) * t221 + t247 * t319;
t220 = t254 * t319 - t282 * t322;
t194 = -qJD(4) * t220 + t247 * t322;
t75 = t193 * pkin(4) - t194 * qJ(5) - t221 * qJD(5) + t161;
t24 = t311 * t75 + t315 * t41;
t424 = Ifges(4,4) * t249;
t422 = Ifges(5,4) * t319;
t421 = Ifges(5,4) * t322;
t420 = Ifges(6,4) * t311;
t419 = Ifges(6,4) * t315;
t417 = Ifges(6,5) * t315;
t416 = Ifges(6,6) * t311;
t408 = t246 * mrSges(4,3);
t173 = -t282 * pkin(3) - t188;
t108 = t220 * pkin(4) - t221 * qJ(5) + t173;
t102 = t319 * t167 + t322 * t174;
t89 = -qJ(5) * t324 + t102;
t50 = t311 * t108 + t315 * t89;
t402 = t144 * t323;
t399 = t246 * t322;
t389 = t316 * (-mrSges(3,2) * t403 + mrSges(3,3) * t392) * qJD(1);
t152 = t294 * t215;
t385 = t152 + t284;
t153 = t337 * t215;
t384 = t153 + t283;
t115 = -mrSges(6,1) * t355 + mrSges(6,2) * t185;
t187 = mrSges(5,1) * t244 - mrSges(5,3) * t216;
t382 = t187 - t115;
t381 = -t493 + t494;
t371 = pkin(10) * t373;
t366 = Ifges(5,5) * t180 - Ifges(5,6) * t181 + Ifges(5,3) * t237;
t17 = -t48 * mrSges(7,1) + t47 * mrSges(7,2);
t357 = t373 / 0.2e1;
t23 = -t311 * t41 + t315 * t75;
t49 = t315 * t108 - t311 * t89;
t101 = t167 * t322 - t319 * t174;
t349 = mrSges(6,1) * t311 + mrSges(6,2) * t315;
t348 = Ifges(5,1) * t322 - t422;
t347 = Ifges(6,1) * t315 - t420;
t346 = -Ifges(5,2) * t319 + t421;
t345 = -Ifges(6,2) * t311 + t419;
t344 = Ifges(5,5) * t322 - Ifges(5,6) * t319;
t343 = -t416 + t417;
t341 = -t21 * t311 + t22 * t315;
t197 = t221 * t315 - t311 * t324;
t28 = pkin(5) * t220 - pkin(11) * t197 + t49;
t196 = -t221 * t311 - t315 * t324;
t34 = pkin(11) * t196 + t50;
t8 = t28 * t321 - t318 * t34;
t9 = t28 * t318 + t321 * t34;
t121 = t196 * t321 - t197 * t318;
t122 = t196 * t318 + t197 * t321;
t338 = -(-qJ(2) * t363 + t303) * t312 + t278 * t316;
t60 = -t319 * t160 - t167 * t374 - t174 * t373 + t199 * t322;
t90 = pkin(4) * t324 - t101;
t334 = t150 * (mrSges(5,1) * t319 + mrSges(5,2) * t322);
t52 = -pkin(4) * t248 - t60;
t285 = (mrSges(3,1) * t403 - mrSges(3,3) * t396) * qJD(1);
t310 = -pkin(5) * t315 - pkin(4);
t296 = t365 * t319;
t280 = t337 * t319;
t279 = t294 * t319;
t272 = -pkin(10) * t397 + t291;
t243 = Ifges(4,4) * t246;
t218 = -mrSges(4,2) * t325 + t408;
t206 = -mrSges(4,1) * t246 + mrSges(4,2) * t249;
t200 = mrSges(4,1) * t237 + mrSges(4,2) * t236;
t192 = Ifges(4,1) * t249 + t325 * Ifges(4,5) + t243;
t191 = Ifges(4,2) * t246 + t325 * Ifges(4,6) + t424;
t186 = -mrSges(5,2) * t244 - mrSges(5,3) * t215;
t155 = t194 * t315 + t248 * t311;
t154 = -t194 * t311 + t248 * t315;
t142 = -mrSges(5,2) * t237 - mrSges(5,3) * t181;
t125 = t216 * Ifges(5,5) - t215 * Ifges(5,6) + t244 * Ifges(5,3);
t124 = mrSges(6,1) * t215 - mrSges(6,3) * t185;
t123 = -mrSges(6,2) * t215 + mrSges(6,3) * t355;
t114 = mrSges(5,1) * t181 + mrSges(5,2) * t180;
t98 = mrSges(6,1) * t181 - mrSges(6,3) * t138;
t97 = -mrSges(6,2) * t181 + mrSges(6,3) * t137;
t93 = t185 * Ifges(6,1) + Ifges(6,4) * t355 + Ifges(6,5) * t215;
t88 = mrSges(7,1) * t214 - mrSges(7,3) * t113;
t87 = -mrSges(7,2) * t214 + mrSges(7,3) * t507;
t77 = -pkin(5) * t401 + t85;
t76 = -pkin(5) * t196 + t90;
t71 = Ifges(6,4) * t138 + Ifges(6,2) * t137 + Ifges(6,6) * t181;
t64 = -pkin(5) * t355 + t79;
t61 = -mrSges(7,1) * t507 + mrSges(7,2) * t113;
t58 = -qJD(6) * t122 + t154 * t321 - t155 * t318;
t57 = qJD(6) * t121 + t154 * t318 + t155 * t321;
t31 = -mrSges(7,2) * t181 + mrSges(7,3) * t48;
t30 = mrSges(7,1) * t181 - mrSges(7,3) * t47;
t29 = -pkin(5) * t154 + t52;
t26 = -pkin(5) * t137 + t37;
t18 = pkin(11) * t154 + t24;
t14 = pkin(5) * t193 - pkin(11) * t155 + t23;
t4 = -qJD(6) * t9 + t14 * t321 - t18 * t318;
t3 = qJD(6) * t8 + t14 * t318 + t18 * t321;
t12 = [(Ifges(6,1) * t155 + Ifges(6,4) * t154) * t442 + (Ifges(6,1) * t197 + Ifges(6,4) * t196) * t447 + (0.2e1 * t389 + m(3) * ((t316 * t380 + (qJ(2) * t396 - t307) * t312) * qJD(1) + t338)) * t378 + (-m(4) * t175 - t467) * t161 + (mrSges(4,2) * t336 + mrSges(4,3) * t144 + Ifges(4,1) * t236 - Ifges(4,4) * t237) * t254 + (Ifges(5,4) * t194 + Ifges(5,6) * t248) * t438 + (Ifges(7,1) * t57 + Ifges(7,4) * t58) * t452 + (Ifges(7,1) * t122 + Ifges(7,4) * t121) * t462 + t194 * t509 - (mrSges(4,1) * t336 + t366 / 0.2e1 + Ifges(5,6) * t445 + Ifges(5,5) * t446 + Ifges(5,3) * t435 - t143 * mrSges(4,3) - Ifges(4,4) * t236 + Ifges(4,2) * t237 + t466) * t324 + (Ifges(6,5) * t197 + Ifges(7,5) * t122 + Ifges(6,6) * t196 + Ifges(7,6) * t121) * t444 + (-t42 * mrSges(5,3) - t518 + t144 * mrSges(5,1) + t496 / 0.2e1 + t513) * t220 + t216 * (Ifges(5,1) * t194 + Ifges(5,5) * t248) / 0.2e1 + (t381 / 0.2e1 - t492 + t494 / 0.2e1 - t493 / 0.2e1 - t144 * mrSges(4,1)) * t282 + (t154 * t36 - t155 * t35 + t196 * t22 - t197 * t21) * mrSges(6,3) + t325 * (Ifges(4,5) * t247 - Ifges(4,6) * t248) / 0.2e1 + t244 * (Ifges(5,5) * t194 + Ifges(5,3) * t248) / 0.2e1 + (-t175 * t247 - t176 * t248 - t188 * t236 - t189 * t237) * mrSges(4,3) + (mrSges(5,2) * t144 - mrSges(5,3) * t43 + 0.2e1 * t456) * t221 + (t1 * t121 - t122 * t2 - t5 * t57 + t58 * t6) * mrSges(7,3) + (Ifges(7,5) * t57 + Ifges(7,6) * t58) * t439 + (Ifges(6,5) * t155 + Ifges(6,6) * t154) * t437 + (Ifges(7,4) * t122 + Ifges(7,2) * t121) * t461 + t248 * t499 - t194 * t511 + m(6) * (t21 * t49 + t22 * t50 + t23 * t35 + t24 * t36 + t37 * t90 + t52 * t79) + m(7) * (t1 * t9 + t2 * t8 + t26 * t76 + t29 * t64 + t3 * t6 + t4 * t5) + m(5) * (t101 * t43 + t102 * t42 + t144 * t173 + t59 * t85 + t60 * t84) + (Ifges(6,5) * t442 - Ifges(5,2) * t438 + Ifges(6,6) * t443 + Ifges(6,3) * t437 + t504 - t512 - t514) * t193 - 0.2e1 * t285 * t362 + t122 * t463 + t121 * t464 + t197 * t458 + t57 * t459 + t58 * t460 + t194 * t449 + (Ifges(4,1) * t247 - Ifges(4,4) * t248) * t432 + t248 * t125 / 0.2e1 + t246 * (Ifges(4,4) * t247 - Ifges(4,2) * t248) / 0.2e1 + t208 * (mrSges(4,1) * t248 + mrSges(4,2) * t247) - t248 * t191 / 0.2e1 + t247 * t192 / 0.2e1 + t217 * t200 + t160 * t218 + t196 * t71 / 0.2e1 + t37 * (-mrSges(6,1) * t196 + mrSges(6,2) * t197) + t59 * t186 + t60 * t187 + t173 * t114 + (Ifges(6,4) * t155 + Ifges(6,2) * t154) * t443 + (Ifges(6,4) * t197 + Ifges(6,2) * t196) * t448 - t248 * t498 + (Ifges(7,4) * t57 + Ifges(7,2) * t58) * t454 + t8 * t30 + t9 * t31 + t206 * t353 + m(4) * (t143 * t189 - t144 * t188 + t160 * t176 + (qJD(1) * t217 + t208) * t353) + t29 * t61 + t64 * (-mrSges(7,1) * t58 + mrSges(7,2) * t57) + t76 * t17 + t3 * t87 + t4 * t88 + t90 * t82 + t50 * t97 + t49 * t98 + t52 * t115 + t26 * (-mrSges(7,1) * t121 + mrSges(7,2) * t122) + t24 * t123 + t23 * t124 + t101 * t141 + t102 * t142 + t154 * t92 / 0.2e1 + t155 * t93 / 0.2e1 + t79 * (-mrSges(6,1) * t154 + mrSges(6,2) * t155); t478 * t87 + t476 * t124 - m(4) * (-t175 * t268 + t176 * t269) + (t17 + t491) * t287 + t477 * t88 + t475 * t123 + t474 * t186 - t383 * t268 + t317 * t200 + t288 * t142 + t259 * t97 - t269 * t218 + t258 * t98 + t211 * t30 + t212 * t31 + (-m(3) * t338 + t312 * t285 - t389) * t379 + t473 * (t61 - t382) + (t1 * t212 + t2 * t211 + t26 * t287 + t473 * t64 + t477 * t5 + t478 * t6) * m(7) + (t21 * t258 + t22 * t259 + t287 * t37 + t35 * t476 + t36 * t475 + t473 * t79) * m(6) + (-t150 * t268 - t43 * t287 + t42 * t288 - t473 * t84 + t474 * t85) * m(5) + (m(5) * (t150 * t377 - t402) - t206 * t363 - t323 * t114 + (-t236 * t323 - t237 * t320) * mrSges(4,3) + (t218 * t323 + t320 * t383) * qJD(3) + (t143 * t320 - t175 * t377 + t176 * t376 - t208 * t363 + t317 * t351 - t402) * m(4)) * t313; (t35 * (mrSges(6,1) * t319 - mrSges(6,3) * t390) + t36 * (-mrSges(6,2) * t319 - mrSges(6,3) * t397) + t334) * qJD(4) + (-mrSges(5,1) * t322 + mrSges(5,2) * t319 - mrSges(4,1)) * t144 + (-t202 / 0.2e1 + t315 * t357) * t93 + (t408 - t218) * t175 + t422 * t445 + t491 * pkin(10) * t319 + (t357 - t399 / 0.2e1) * t127 + (t407 + t467) * t176 + (Ifges(6,5) * t202 + Ifges(6,6) * t201) * t438 + (t371 - t107) * t115 - (Ifges(4,1) * t246 + t125 - t424) * t249 / 0.2e1 - (-Ifges(4,2) * t249 + t192 + t243) * t246 / 0.2e1 - t355 * (Ifges(6,4) * t202 + Ifges(6,2) * t201) / 0.2e1 + (t185 * (Ifges(6,5) * t319 + t322 * t347) + t355 * (Ifges(6,6) * t319 + t322 * t345) + t244 * t344 + t216 * t348) * qJD(4) / 0.2e1 + (t427 - t495 / 0.2e1 + Ifges(7,5) * t453 + Ifges(7,6) * t455 + Ifges(6,3) * t438 + Ifges(7,3) * t440 + t508) * t400 - t325 * (Ifges(4,5) * t246 - Ifges(4,6) * t249) / 0.2e1 + (-t201 * t36 + t202 * t35 - t21 * t391 - t22 * t398) * mrSges(6,3) + t485 * t123 + t486 * t124 + (pkin(10) * t472 - t107 * t79 + t21 * t272 + t22 * t273 + t35 * t486 + t36 * t485) * m(6) + t489 * t87 + (t1 * t210 + t2 * t209 + t26 * t296 + t484 * t64 + t489 * t6 + t490 * t5) * m(7) + t490 * t88 + (Ifges(7,4) * t230 + Ifges(7,2) * t231) * t454 + (Ifges(7,4) * t129 + Ifges(7,2) * t128) * t455 + (Ifges(6,3) * t502 + t343 * t444 + t345 * t448 + t347 * t447 + t456 + t505) * t319 + (t399 * t84 + t471) * mrSges(5,3) - t185 * (Ifges(6,1) * t202 + Ifges(6,4) * t201) / 0.2e1 + (-t371 - t118) * t187 + (-Ifges(7,4) * t280 - Ifges(7,2) * t279) * t461 + (-Ifges(7,5) * t280 - Ifges(7,6) * t279) * t444 + (-Ifges(7,1) * t280 - Ifges(7,4) * t279) * t462 + t26 * (mrSges(7,1) * t279 - mrSges(7,2) * t280) + (Ifges(7,5) * t230 + Ifges(7,6) * t231) * t439 + (Ifges(7,5) * t129 + Ifges(7,6) * t128) * t440 + t249 * t498 + t470 * t373 + (-t118 * t84 - t119 * t85 - pkin(3) * t144 + ((-t319 * t85 - t322 * t84) * qJD(4) + t471) * pkin(10)) * m(5) + t472 * t349 + (mrSges(7,1) * t480 - mrSges(7,2) * t479) * t64 + (-t1 * t279 + t2 * t280 + t479 * t5 - t480 * t6) * mrSges(7,3) + t381 + (Ifges(7,1) * t230 + Ifges(7,4) * t231) * t452 + (Ifges(7,1) * t129 + Ifges(7,4) * t128) * t453 + (-t372 - t119) * t186 + (pkin(10) * t142 + t343 * t502 - t513) * t322 - t71 * t398 / 0.2e1 - t346 * t375 / 0.2e1 - t280 * t463 - t279 * t464 + t201 * t457 + t391 * t458 + t230 * t459 + t231 * t460 + t191 * t432 + (Ifges(5,6) * t249 + t246 * t346) * t437 + t296 * t17 + t272 * t98 + t273 * t97 - t208 * (mrSges(4,1) * t249 + mrSges(4,2) * t246) + t209 * t30 + t210 * t31 - t79 * (-mrSges(6,1) * t201 + mrSges(6,2) * t202) - t496 * t322 / 0.2e1 - t249 * t499 - t492 + (-t427 - t468 + t504) * t374 + t421 * t446 + t484 * t61 - t246 * t334 - pkin(3) * t114 - t128 * t54 / 0.2e1 - t129 * t55 / 0.2e1 - t244 * (Ifges(5,3) * t249 + t246 * t344) / 0.2e1 - t216 * (Ifges(5,5) * t249 + t246 * t348) / 0.2e1; (-t35 * t68 - t36 * t69 - t79 * t85 - pkin(4) * t37 + (-t311 * t35 + t315 * t36) * qJD(5) + t341 * qJ(5)) * m(6) + (-t411 / 0.2e1 - t53 / 0.2e1 - t91 / 0.2e1 - t415 / 0.2e1 - t414 / 0.2e1 + t508 + t512) * t216 + (Ifges(6,5) * t311 + Ifges(7,5) * t294 + Ifges(6,6) * t315 - Ifges(7,6) * t337) * t444 + (-t1 * t337 - t2 * t294 + t384 * t5 - t385 * t6) * mrSges(7,3) + (Ifges(7,4) * t294 - Ifges(7,2) * t337) * t461 + (Ifges(7,1) * t294 - Ifges(7,4) * t337) * t462 + t26 * (mrSges(7,1) * t337 + mrSges(7,2) * t294) - t337 * t464 + t487 * t88 + (t1 * t266 + t2 * t265 + t26 * t310 + t487 * t5 + t488 * t6 - t64 * t77) * m(7) + t488 * t87 + (-t311 * t98 + t315 * t97) * qJ(5) + (t79 * t349 + t345 * t443 + t347 * t442 + t410 / 0.2e1 - t213 / 0.2e1 + t509 + t449 + t93 * t431 + (t417 / 0.2e1 - t416 / 0.2e1) * t215 + (-t311 * t36 - t315 * t35) * mrSges(6,3) + (-Ifges(6,3) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t216 + t470) * t215 + (-t283 / 0.2e1 - t153 / 0.2e1) * t55 + (-Ifges(7,1) * t283 - Ifges(7,4) * t284) * t452 + (-Ifges(7,4) * t283 - Ifges(7,2) * t284) * t454 + (-Ifges(7,5) * t283 - Ifges(7,6) * t284) * t439 + (-t284 / 0.2e1 - t152 / 0.2e1) * t54 + (t123 * t315 - t124 * t311) * qJD(5) + t466 + t366 + (mrSges(7,1) * t385 - mrSges(7,2) * t384) * t64 + t382 * t85 + t294 * t463 + (Ifges(7,1) * t153 + Ifges(7,4) * t152) * t453 + (Ifges(7,4) * t153 + Ifges(7,2) * t152) * t455 + t311 * t458 + (Ifges(7,5) * t153 + Ifges(7,6) * t152) * t440 + (Ifges(6,1) * t311 + t419) * t447 + (Ifges(6,2) * t315 + t420) * t448 + t71 * t431 + t37 * (-mrSges(6,1) * t315 + mrSges(6,2) * t311) + t310 * t17 + t265 * t30 + t266 * t31 - t84 * t186 - t77 * t61 - pkin(4) * t82 - t69 * t123 - t68 * t124 + t341 * mrSges(6,3); t113 * t88 - t507 * t87 - t355 * t123 + t185 * t124 + t17 + t82 + (t113 * t5 - t507 * t6 + t26) * m(7) + (t185 * t35 - t355 * t36 + t37) * m(6); -t64 * (mrSges(7,1) * t113 + mrSges(7,2) * t507) + (Ifges(7,1) * t507 - t418) * t453 + t54 * t452 + (Ifges(7,5) * t507 - Ifges(7,6) * t113) * t440 - t5 * t87 + t6 * t88 + (t113 * t6 + t5 * t507) * mrSges(7,3) + t11 + (-Ifges(7,2) * t113 + t109 + t55) * t455 + t469;];
tauc  = t12(:);
