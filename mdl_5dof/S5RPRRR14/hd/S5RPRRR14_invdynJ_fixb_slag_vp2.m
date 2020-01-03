% Calculate vector of inverse dynamics joint torques for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR14_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:16:54
% EndTime: 2019-12-31 19:17:36
% DurationCPUTime: 21.25s
% Computational Cost: add. (16135->709), mult. (51743->1013), div. (0->0), fcn. (44289->14), ass. (0->342)
t256 = sin(pkin(11));
t258 = sin(pkin(5));
t259 = cos(pkin(11));
t267 = cos(qJ(3));
t260 = cos(pkin(6));
t264 = sin(qJ(3));
t372 = t260 * t264;
t280 = t258 * (-t256 * t372 + t259 * t267);
t206 = qJD(1) * t280;
t257 = sin(pkin(6));
t361 = qJD(3) * t267;
t338 = t257 * t361;
t497 = -t206 + t338;
t261 = cos(pkin(5));
t371 = t260 * t267;
t293 = -t256 * t264 + t259 * t371;
t377 = t257 * t267;
t275 = t258 * t293 + t261 * t377;
t180 = t275 * qJD(1);
t175 = -t180 + qJD(4);
t292 = t256 * t267 + t259 * t372;
t354 = qJD(1) * qJD(3);
t379 = t257 * t261;
t130 = (qJDD(1) * t264 + t267 * t354) * t379 + (qJDD(1) * t292 + t293 * t354) * t258;
t378 = t257 * t264;
t191 = t258 * t292 + t261 * t378;
t183 = t191 * qJD(1);
t376 = t258 * t259;
t222 = -t257 * t376 + t260 * t261;
t213 = qJD(1) * t222 + qJD(3);
t263 = sin(qJ(4));
t266 = cos(qJ(4));
t140 = -t183 * t263 + t213 * t266;
t212 = qJDD(1) * t222 + qJDD(3);
t75 = qJD(4) * t140 + t130 * t266 + t212 * t263;
t430 = t75 / 0.2e1;
t496 = Ifges(5,4) * t430;
t229 = -t260 * t266 + t263 * t378;
t381 = t256 * t258;
t341 = qJD(1) * t381;
t327 = t257 * t341;
t480 = -qJD(4) * t229 - t263 * t327 + t266 * t497;
t279 = (t256 * t371 + t259 * t264) * t258;
t205 = qJD(1) * t279;
t362 = qJD(3) * t264;
t339 = t257 * t362;
t495 = -t205 + t339;
t413 = sin(qJ(1));
t253 = t413 * t256;
t268 = cos(qJ(1));
t365 = t268 * t259;
t321 = t261 * t365 - t253;
t374 = t258 * t268;
t494 = -t257 * t374 + t260 * t321;
t342 = t413 * t259;
t366 = t268 * t256;
t223 = t261 * t366 + t342;
t152 = t223 * t267 + t264 * t494;
t193 = t257 * t321 + t260 * t374;
t114 = t152 * t266 - t193 * t263;
t493 = -t152 * t263 - t193 * t266;
t141 = t183 * t266 + t213 * t263;
t76 = -qJD(4) * t141 - t130 * t263 + t212 * t266;
t429 = t76 / 0.2e1;
t131 = t258 * (qJDD(1) * t293 - t292 * t354) - (-qJDD(1) * t267 + t264 * t354) * t379;
t128 = qJDD(4) - t131;
t423 = t128 / 0.2e1;
t262 = sin(qJ(5));
t265 = cos(qJ(5));
t412 = pkin(1) * t261;
t352 = qJD(1) * t412;
t246 = t259 * t352;
t411 = pkin(2) * t261;
t278 = t411 + (-pkin(8) * t260 - qJ(2)) * t381;
t176 = qJD(1) * t278 + t246;
t382 = t256 * t257;
t208 = (-pkin(2) * t259 - pkin(8) * t382 - pkin(1)) * t258;
t199 = qJD(1) * t208 + qJD(2);
t137 = -t176 * t257 + t199 * t260;
t81 = -pkin(3) * t180 - pkin(9) * t183 + t137;
t248 = qJ(2) * t376;
t218 = qJD(1) * t248 + t256 * t352;
t375 = t258 * t260;
t283 = (t259 * t375 + t379) * pkin(8);
t170 = qJD(1) * t283 + t218;
t385 = t176 * t260;
t300 = t199 * t257 + t385;
t97 = t170 * t267 + t264 * t300;
t83 = pkin(9) * t213 + t97;
t39 = t263 * t81 + t266 * t83;
t36 = pkin(10) * t175 + t39;
t96 = -t170 * t264 + t267 * t300;
t82 = -pkin(3) * t213 - t96;
t48 = -pkin(4) * t140 - pkin(10) * t141 + t82;
t14 = -t262 * t36 + t265 * t48;
t15 = t262 * t48 + t265 * t36;
t446 = t14 * mrSges(6,1) - t15 * mrSges(6,2);
t100 = t141 * t265 + t175 * t262;
t139 = qJD(5) - t140;
t99 = -t141 * t262 + t175 * t265;
t43 = Ifges(6,5) * t100 + Ifges(6,6) * t99 + Ifges(6,3) * t139;
t396 = Ifges(5,4) * t141;
t485 = t175 * Ifges(5,6);
t486 = t140 * Ifges(5,2);
t73 = t396 + t485 + t486;
t476 = -t73 / 0.2e1 + t43 / 0.2e1;
t492 = t446 + t476;
t71 = qJDD(5) - t76;
t432 = t71 / 0.2e1;
t33 = -qJD(5) * t100 + t128 * t265 - t262 * t75;
t438 = t33 / 0.2e1;
t32 = qJD(5) * t99 + t128 * t262 + t265 * t75;
t439 = t32 / 0.2e1;
t231 = (qJ(2) * qJDD(1) + qJD(1) * qJD(2)) * t258;
t348 = qJDD(1) * t412;
t202 = -t231 * t256 + t259 * t348;
t163 = (-pkin(8) * t256 * t375 + t411) * qJDD(1) + t202;
t195 = qJDD(1) * t208 + qJDD(2);
t203 = t231 * t259 + t256 * t348;
t481 = qJD(3) * t385 + qJDD(1) * t283 + t203;
t50 = t267 * (t163 * t260 + t195 * t257) - t170 * t361 - t199 * t339 - t481 * t264;
t47 = -pkin(3) * t212 - t50;
t13 = -pkin(4) * t76 - pkin(10) * t75 + t47;
t359 = qJD(4) * t266;
t360 = qJD(4) * t263;
t49 = t163 * t372 - t170 * t362 + t195 * t378 + t199 * t338 + t267 * t481;
t46 = pkin(9) * t212 + t49;
t121 = -t163 * t257 + t195 * t260;
t61 = -pkin(3) * t131 - pkin(9) * t130 + t121;
t10 = t263 * t61 + t266 * t46 + t359 * t81 - t360 * t83;
t5 = pkin(10) * t128 + t10;
t1 = qJD(5) * t14 + t13 * t262 + t265 * t5;
t2 = -qJD(5) * t15 + t13 * t265 - t262 * t5;
t447 = -t2 * mrSges(6,1) + t1 * mrSges(6,2);
t7 = Ifges(6,5) * t32 + Ifges(6,6) * t33 + Ifges(6,3) * t71;
t491 = t447 + 0.2e1 * Ifges(5,2) * t429 + 0.2e1 * Ifges(5,6) * t423 + t496 - t7 / 0.2e1 - Ifges(6,5) * t439 - Ifges(6,6) * t438 - Ifges(6,3) * t432;
t469 = m(6) + m(5);
t350 = t469 + m(3) + m(4);
t132 = pkin(3) * t183 - pkin(9) * t180;
t64 = t132 * t263 + t266 * t96;
t490 = pkin(9) * t360 + pkin(10) * t183 + t64;
t489 = -pkin(9) * qJD(5) * t266 - t97 + t175 * (pkin(4) * t263 - pkin(10) * t266);
t488 = t82 * mrSges(5,1);
t487 = t82 * mrSges(5,2);
t398 = mrSges(5,3) * t141;
t102 = mrSges(5,1) * t175 - t398;
t62 = -mrSges(6,1) * t99 + mrSges(6,2) * t100;
t484 = t62 - t102;
t230 = t260 * t263 + t266 * t378;
t294 = -t230 * t265 + t262 * t377;
t483 = qJD(5) * t294 - t262 * t480 + t265 * t495;
t200 = -t230 * t262 - t265 * t377;
t482 = qJD(5) * t200 + t262 * t495 + t265 * t480;
t226 = t256 * t412 + t248;
t184 = t283 + t226;
t252 = t259 * t412;
t192 = t252 + t278;
t299 = t192 * t260 + t208 * t257;
t103 = -t184 * t264 + t267 * t299;
t479 = qJD(4) * t230 + t263 * t497 + t266 * t327;
t284 = t261 * t342 + t366;
t343 = t258 * t413;
t478 = t257 * t284 + t260 * t343;
t369 = t262 * t266;
t122 = -t180 * t369 + t183 * t265;
t477 = t262 * t359 + t122;
t317 = -mrSges(6,1) * t265 + mrSges(6,2) * t262;
t287 = m(6) * pkin(4) - t317;
t319 = -mrSges(5,1) * t266 + mrSges(5,2) * t263;
t349 = m(6) * pkin(10) + mrSges(6,3);
t452 = t263 * t349 + t266 * t287 + mrSges(4,1) - t319;
t475 = pkin(3) * t469 + t452;
t474 = -t223 * t264 + t267 * t494;
t473 = Ifges(5,1) * t430 + Ifges(5,5) * t423;
t440 = Ifges(5,4) * t429 + t473;
t421 = t139 / 0.2e1;
t424 = t100 / 0.2e1;
t427 = t99 / 0.2e1;
t472 = Ifges(6,5) * t424 + Ifges(6,6) * t427 + Ifges(6,3) * t421;
t38 = -t263 * t83 + t266 * t81;
t468 = t38 * mrSges(5,1);
t467 = t39 * mrSges(5,2);
t12 = -mrSges(6,1) * t33 + mrSges(6,2) * t32;
t53 = mrSges(5,1) * t128 - mrSges(5,3) * t75;
t466 = t12 - t53;
t465 = mrSges(5,1) + t287;
t316 = mrSges(6,1) * t262 + mrSges(6,2) * t265;
t464 = -m(6) * pkin(9) + mrSges(4,2) - t316;
t331 = mrSges(5,2) - t349;
t239 = -pkin(4) * t266 - pkin(10) * t263 - pkin(3);
t358 = qJD(5) * t262;
t463 = -t239 * t358 + t262 * t490 + t265 * t489;
t356 = qJD(5) * t265;
t462 = t239 * t356 + t262 * t489 - t265 * t490;
t400 = mrSges(4,3) * t183;
t388 = mrSges(4,1) * t213 + mrSges(5,1) * t140 - mrSges(5,2) * t141 - t400;
t458 = t263 * t356 + t477;
t367 = t265 * t266;
t123 = t180 * t367 + t183 * t262;
t357 = qJD(5) * t263;
t457 = t262 * t357 - t265 * t359 + t123;
t456 = -t257 * t343 + t260 * t284;
t11 = -qJD(4) * t39 - t263 * t46 + t266 * t61;
t455 = t10 * t266 - t11 * t263;
t454 = t1 * t265 - t2 * t262;
t451 = t256 * ((-mrSges(4,1) * t180 + mrSges(4,2) * t183) * t257 - (mrSges(3,1) * t261 - mrSges(3,3) * t381) * qJD(1)) + t259 * (-mrSges(3,2) * t261 + mrSges(3,3) * t376) * qJD(1);
t450 = t11 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,5) * t75 + Ifges(5,6) * t76 + Ifges(5,3) * t128;
t181 = t275 * qJD(3);
t182 = t191 * qJD(3);
t326 = qJD(2) * t257 * t381;
t120 = pkin(3) * t182 - pkin(9) * t181 + t326;
t142 = -t192 * t257 + t208 * t260;
t188 = t275 * pkin(3);
t334 = pkin(9) * t191 + t188;
t91 = t142 - t334;
t173 = t267 * t184;
t104 = t192 * t372 + t208 * t378 + t173;
t95 = pkin(9) * t222 + t104;
t404 = t263 * t91 + t266 * t95;
t85 = qJD(2) * t280 + qJD(3) * t103;
t23 = -qJD(4) * t404 + t120 * t266 - t263 * t85;
t449 = -m(5) * pkin(9) + t464;
t448 = t50 * mrSges(4,1) - t49 * mrSges(4,2) + Ifges(4,5) * t130 + Ifges(4,6) * t131 + Ifges(4,3) * t212;
t422 = -t139 / 0.2e1;
t425 = -t100 / 0.2e1;
t428 = -t99 / 0.2e1;
t444 = Ifges(6,5) * t425 + Ifges(6,6) * t428 + Ifges(6,3) * t422 - t446;
t443 = 0.2e1 * t261;
t8 = Ifges(6,4) * t32 + Ifges(6,2) * t33 + Ifges(6,6) * t71;
t442 = t8 / 0.2e1;
t441 = Ifges(6,1) * t439 + Ifges(6,4) * t438 + Ifges(6,5) * t432;
t393 = Ifges(6,4) * t100;
t44 = Ifges(6,2) * t99 + Ifges(6,6) * t139 + t393;
t436 = -t44 / 0.2e1;
t435 = t44 / 0.2e1;
t98 = Ifges(6,4) * t99;
t45 = Ifges(6,1) * t100 + Ifges(6,5) * t139 + t98;
t434 = -t45 / 0.2e1;
t433 = t45 / 0.2e1;
t420 = -t140 / 0.2e1;
t419 = -t141 / 0.2e1;
t418 = t141 / 0.2e1;
t417 = -t175 / 0.2e1;
t414 = t183 / 0.2e1;
t6 = -pkin(4) * t128 - t11;
t407 = t263 * t6;
t403 = mrSges(3,1) * t259;
t402 = mrSges(3,2) * t256;
t401 = mrSges(4,3) * t180;
t399 = mrSges(5,3) * t140;
t397 = Ifges(4,4) * t183;
t395 = Ifges(5,4) * t263;
t394 = Ifges(5,4) * t266;
t392 = Ifges(6,4) * t262;
t391 = Ifges(6,4) * t265;
t387 = t140 * t262;
t386 = t140 * t265;
t384 = t180 * t263;
t383 = t180 * t266;
t370 = t262 * t263;
t368 = t263 * t265;
t363 = pkin(1) * t268 + qJ(2) * t343;
t353 = qJDD(1) * t258;
t351 = pkin(9) * t359;
t336 = t359 / 0.2e1;
t335 = -t357 / 0.2e1;
t325 = -pkin(1) * t413 + qJ(2) * t374;
t149 = t191 * t263 - t222 * t266;
t150 = t191 * t266 + t222 * t263;
t320 = mrSges(5,1) * t149 + mrSges(5,2) * t150;
t111 = -t150 * t262 - t265 * t275;
t112 = t150 * t265 - t262 * t275;
t318 = mrSges(6,1) * t111 - mrSges(6,2) * t112;
t315 = Ifges(5,1) * t266 - t395;
t314 = Ifges(6,1) * t265 - t392;
t313 = Ifges(6,1) * t262 + t391;
t312 = -Ifges(5,2) * t263 + t394;
t311 = -Ifges(6,2) * t262 + t391;
t310 = Ifges(6,2) * t265 + t392;
t309 = Ifges(5,5) * t266 - Ifges(5,6) * t263;
t308 = Ifges(6,5) * t265 - Ifges(6,6) * t262;
t307 = Ifges(6,5) * t262 + Ifges(6,6) * t265;
t41 = -pkin(10) * t275 + t404;
t94 = -pkin(3) * t222 - t103;
t58 = pkin(4) * t149 - pkin(10) * t150 + t94;
t19 = t262 * t58 + t265 * t41;
t18 = -t262 * t41 + t265 * t58;
t51 = -t263 * t95 + t266 * t91;
t304 = -pkin(9) * t469 + mrSges(4,2) - mrSges(5,3);
t63 = t132 * t266 - t263 * t96;
t298 = -(-qJ(2) * t341 + t246) * t256 + t218 * t259;
t35 = -pkin(4) * t175 - t38;
t291 = t35 * t316;
t22 = t120 * t263 + t266 * t85 + t359 * t91 - t360 * t95;
t274 = -pkin(2) * t223 + pkin(8) * t193 + t325;
t224 = -t253 * t261 + t365;
t271 = pkin(2) * t224 + pkin(8) * t478 + t363;
t156 = t224 * t267 - t264 * t456;
t269 = pkin(3) * t156 + t271;
t86 = qJD(2) * t279 + (t264 * t299 + t173) * qJD(3);
t247 = -pkin(1) * t353 + qJDD(2);
t240 = t353 * t402;
t225 = -qJ(2) * t381 + t252;
t215 = pkin(9) * t367 + t239 * t262;
t214 = -pkin(9) * t369 + t239 * t265;
t174 = Ifges(4,4) * t180;
t155 = t224 * t264 + t267 * t456;
t143 = -mrSges(4,2) * t213 + t401;
t138 = Ifges(5,4) * t140;
t118 = t156 * t266 + t263 * t478;
t117 = t156 * t263 - t266 * t478;
t110 = -qJD(4) * t149 + t181 * t266;
t109 = qJD(4) * t150 + t181 * t263;
t108 = t183 * Ifges(4,1) + t213 * Ifges(4,5) + t174;
t107 = t180 * Ifges(4,2) + t213 * Ifges(4,6) + t397;
t106 = -mrSges(4,2) * t212 + mrSges(4,3) * t131;
t105 = mrSges(4,1) * t212 - mrSges(4,3) * t130;
t101 = -mrSges(5,2) * t175 + t399;
t89 = pkin(4) * t141 - pkin(10) * t140;
t79 = t118 * t265 + t155 * t262;
t78 = -t118 * t262 + t155 * t265;
t77 = -mrSges(4,1) * t131 + mrSges(4,2) * t130;
t74 = t141 * Ifges(5,1) + t175 * Ifges(5,5) + t138;
t72 = Ifges(5,5) * t141 + t140 * Ifges(5,6) + t175 * Ifges(5,3);
t67 = mrSges(6,1) * t139 - mrSges(6,3) * t100;
t66 = -mrSges(6,2) * t139 + mrSges(6,3) * t99;
t60 = qJD(5) * t111 + t110 * t265 + t182 * t262;
t59 = -qJD(5) * t112 - t110 * t262 + t182 * t265;
t55 = -pkin(4) * t183 - t63;
t54 = -mrSges(5,2) * t128 + mrSges(5,3) * t76;
t40 = pkin(4) * t275 - t51;
t37 = pkin(4) * t109 - pkin(10) * t110 + t86;
t34 = -mrSges(5,1) * t76 + mrSges(5,2) * t75;
t27 = t262 * t89 + t265 * t38;
t26 = -t262 * t38 + t265 * t89;
t21 = -pkin(4) * t182 - t23;
t20 = pkin(10) * t182 + t22;
t17 = -mrSges(6,2) * t71 + mrSges(6,3) * t33;
t16 = mrSges(6,1) * t71 - mrSges(6,3) * t32;
t4 = -qJD(5) * t19 - t20 * t262 + t265 * t37;
t3 = qJD(5) * t18 + t20 * t265 + t262 * t37;
t9 = [-t6 * t318 + t47 * t320 + m(3) * (t202 * t225 + t203 * t226 + (-pkin(1) * t247 + qJD(2) * t298) * t258) + m(5) * (t10 * t404 + t11 * t51 + t22 * t39 + t23 * t38 + t47 * t94 + t82 * t86) + t404 * t54 + (t1 * t111 - t112 * t2 - t14 * t60 + t15 * t59) * mrSges(6,3) + (-t39 * mrSges(5,3) - t486 / 0.2e1 - t485 / 0.2e1 - Ifges(5,4) * t418 + t488 + t472 + t492) * t109 - t182 * t467 + m(4) * (t103 * t50 + t104 * t49 + t121 * t142 + t137 * t326 + t85 * t97 - t86 * t96) + (-m(3) * t325 + t223 * mrSges(3,1) + t321 * mrSges(3,2) - mrSges(3,3) * t374 + t413 * mrSges(2,1) + t268 * mrSges(2,2) - m(4) * t274 + t152 * mrSges(4,1) - t193 * mrSges(4,3) + t465 * t114 + (t304 - t316) * t474 + t331 * t493 - t469 * (-pkin(3) * t152 + t274)) * g(1) + (-pkin(1) * t240 + t247 * (t402 - t403) + (-t202 * t256 + t203 * t259) * mrSges(3,3) + t451 * qJD(2)) * t258 + t448 * t222 + (Ifges(6,4) * t60 + Ifges(6,2) * t59) * t427 + (Ifges(6,4) * t112 + Ifges(6,2) * t111) * t438 + (-mrSges(5,3) * t10 - t491 - t496) * t149 + t140 * (Ifges(5,4) * t110 + Ifges(5,6) * t182) / 0.2e1 + (Ifges(6,1) * t60 + Ifges(6,4) * t59) * t424 + (Ifges(6,1) * t112 + Ifges(6,4) * t111) * t439 - t38 * mrSges(5,3) * t110 + t175 * (Ifges(5,5) * t110 + Ifges(5,3) * t182) / 0.2e1 + (-t181 * t96 - t182 * t97) * mrSges(4,3) + t110 * t487 + (t202 * mrSges(3,1) - t203 * mrSges(3,2)) * t261 + (Ifges(6,5) * t60 + Ifges(6,6) * t59) * t421 + (Ifges(6,5) * t112 + Ifges(6,6) * t111) * t432 - (mrSges(4,1) * t121 - mrSges(4,3) * t49 - Ifges(4,4) * t130 - Ifges(4,2) * t131 - Ifges(4,6) * t212 + t450) * t275 + (-m(4) * t271 - t156 * mrSges(4,1) - t478 * mrSges(4,3) - m(3) * t363 - t224 * mrSges(3,1) + t284 * mrSges(3,2) - mrSges(3,3) * t343 - t268 * mrSges(2,1) + t413 * mrSges(2,2) - m(6) * (t118 * pkin(4) + t269) - t79 * mrSges(6,1) - t78 * mrSges(6,2) - m(5) * t269 - t118 * mrSges(5,1) + t304 * t155 + t331 * t117) * g(2) + m(6) * (t1 * t19 + t14 * t4 + t15 * t3 + t18 * t2 + t21 * t35 + t40 * t6) + (-mrSges(5,3) * t11 + 0.2e1 * t440) * t150 + (Ifges(5,1) * t110 + Ifges(5,5) * t182) * t418 + (mrSges(4,2) * t121 - mrSges(4,3) * t50 + Ifges(4,1) * t130 + Ifges(4,4) * t131 + Ifges(4,5) * t212) * t191 + (Ifges(2,3) + (mrSges(3,1) * t225 - mrSges(3,2) * t226 + Ifges(3,3) * t261) * t261 + ((-mrSges(3,3) * t225 + Ifges(3,1) * t381 + Ifges(3,5) * t443) * t256 + (t226 * mrSges(3,3) + Ifges(3,6) * t443 + (mrSges(3,1) * pkin(1) + 0.2e1 * Ifges(3,4) * t256 + Ifges(3,2) * t259) * t258) * t259) * t258) * qJDD(1) + t213 * (Ifges(4,5) * t181 - Ifges(4,6) * t182) / 0.2e1 + t182 * t72 / 0.2e1 - t182 * t107 / 0.2e1 + t137 * (mrSges(4,1) * t182 + mrSges(4,2) * t181) + t180 * (Ifges(4,4) * t181 - Ifges(4,2) * t182) / 0.2e1 + t182 * t468 + (Ifges(4,1) * t181 - Ifges(4,4) * t182) * t414 + t60 * t433 + t59 * t435 + t112 * t441 + t111 * t442 + t18 * t16 + t19 * t17 + t40 * t12 + t51 * t53 + t35 * (-mrSges(6,1) * t59 + mrSges(6,2) * t60) + t21 * t62 + t3 * t66 + t4 * t67 + t94 * t34 + t22 * t101 + t23 * t102 + t103 * t105 + t104 * t106 + t110 * t74 / 0.2e1 - t388 * t86 + t142 * t77 + t85 * t143 + t181 * t108 / 0.2e1; -t206 * t143 + t200 * t16 - t294 * t17 + t230 * t54 + t240 + t483 * t67 + t482 * t66 + t466 * t229 + t388 * t205 + t480 * t101 - t350 * t261 * g(3) - m(4) * (-t205 * t96 + t206 * t97) + m(3) * t247 + (-qJDD(1) * t403 + (-m(4) * t137 * t382 - m(3) * t298 - t451) * qJD(1) + (-g(1) * t413 + g(2) * t268) * t350) * t258 + t484 * t479 + (m(4) * t121 + t77) * t260 + (-t1 * t294 + t14 * t483 + t15 * t482 + t2 * t200 + t229 * t6 + t35 * t479) * m(6) + (t10 * t230 - t11 * t229 - t205 * t82 - t479 * t38 + t480 * t39) * m(5) + (t264 * t106 + (t105 - t34) * t267 + (t143 * t267 - t264 * t388) * qJD(3) + m(5) * (-t267 * t47 + t362 * t82) + m(4) * (t264 * t49 + t267 * t50 + t361 * t97 - t362 * t96)) * t257; (t335 * t44 + t336 * t45) * t265 + t47 * t319 + (t140 * t312 + t141 * t315 + t175 * t309) * qJD(4) / 0.2e1 - (Ifges(4,1) * t180 - t397 + t72) * t183 / 0.2e1 - (-Ifges(4,2) * t183 + t108 + t174) * t180 / 0.2e1 + (-t307 * t421 - t310 * t427 - t313 * t424) * t357 + t448 + (t444 - t476) * t384 + t466 * pkin(9) * t263 + (pkin(9) * t54 + (t308 * t421 + t311 * t427 + t314 * t424) * qJD(4) + t491) * t266 + (-pkin(9) * t101 + t492) * t360 - t183 * t468 + (-t351 - t63) * t102 + (t155 * t475 + t156 * t449) * g(1) + t477 * t436 + (t458 * mrSges(6,1) - t457 * mrSges(6,2)) * t35 + (-t1 * t370 + t14 * t457 - t15 * t458 - t2 * t368) * mrSges(6,3) + (-g(1) * t156 - g(2) * t152 - g(3) * t191 + (-t360 + t384) * t39 + (-t359 + t383) * t38 + t455) * mrSges(5,3) + (-t38 * t63 - t39 * t64 - pkin(3) * t47 + ((-t263 * t39 - t266 * t38) * qJD(4) + t455) * pkin(9)) * m(5) + (t401 - t143) * t96 + (t336 - t383 / 0.2e1) * t74 + (qJD(4) * t472 + t308 * t432 + t311 * t438 + t314 * t439 + t440 + t473) * t263 + (-m(5) * t82 + t388 + t400) * t97 + t462 * t66 + (-t35 * t55 + t1 * t215 + t2 * t214 + (t35 * t359 + t407) * pkin(9) + t462 * t15 + t463 * t14) * m(6) + t463 * t67 + (Ifges(6,5) * t123 + Ifges(6,6) * t122) * t422 + t395 * t429 + t262 * t45 * t335 + (Ifges(6,4) * t123 + Ifges(6,2) * t122) * t428 + t394 * t430 + (t152 * t449 - t474 * t475) * g(2) + (-m(5) * t334 - m(6) * t188 + t191 * t464 - t275 * t452) * g(3) + (Ifges(6,1) * t123 + Ifges(6,4) * t122) * t425 + t175 * t82 * (mrSges(5,1) * t263 + mrSges(5,2) * t266) + (t351 - t55) * t62 - t213 * (Ifges(4,5) * t180 - Ifges(4,6) * t183) / 0.2e1 + t214 * t16 + t215 * t17 - t137 * (mrSges(4,1) * t183 + mrSges(4,2) * t180) + t183 * t467 + t316 * t407 + t107 * t414 + (Ifges(5,3) * t183 + t180 * t309) * t417 + (Ifges(5,5) * t183 + t180 * t315) * t419 + (Ifges(5,6) * t183 + t180 * t312) * t420 + t123 * t434 + t368 * t441 - pkin(3) * t34 - t8 * t370 / 0.2e1 - t64 * t101; t6 * t317 + (t100 * t314 + t139 * t308 + t311 * t99) * qJD(5) / 0.2e1 + (t138 + t74) * t420 + (t331 * t114 - t465 * t493) * g(2) + (-t396 + t43) * t419 + (Ifges(5,1) * t419 + Ifges(5,5) * t417 + t308 * t422 + t311 * t428 + t314 * t425 - t291 - t487) * t140 + (-Ifges(5,2) * t420 - Ifges(5,6) * t417 + t444 - t488) * t141 + (-m(6) * t35 + t398 - t484) * t39 + ((-t358 + t387) * t15 + (-t356 + t386) * t14 + t454) * mrSges(6,3) + (-t66 * t358 - t67 * t356 + m(6) * ((-t14 * t265 - t15 * t262) * qJD(5) + t454) - t262 * t16 + t265 * t17) * pkin(10) + (-pkin(4) * t6 - t14 * t26 - t15 * t27) * m(6) + (t117 * t465 + t118 * t331) * g(1) + (t149 * t287 - t150 * t349 + t320) * g(3) + t450 + qJD(5) * t291 + t73 * t418 + t307 * t432 + t356 * t433 + t386 * t434 + t387 * t435 + t358 * t436 + t310 * t438 + t313 * t439 + t262 * t441 + t265 * t442 - pkin(4) * t12 + (t399 - t101) * t38 - t27 * t66 - t26 * t67; -t35 * (mrSges(6,1) * t100 + mrSges(6,2) * t99) + (Ifges(6,1) * t99 - t393) * t425 + t44 * t424 + (Ifges(6,5) * t99 - Ifges(6,6) * t100) * t422 - t14 * t66 + t15 * t67 - g(1) * (mrSges(6,1) * t78 - mrSges(6,2) * t79) - g(2) * ((-t114 * t262 - t265 * t474) * mrSges(6,1) + (-t114 * t265 + t262 * t474) * mrSges(6,2)) - g(3) * t318 + (t100 * t15 + t14 * t99) * mrSges(6,3) + t7 + (-Ifges(6,2) * t100 + t45 + t98) * t428 - t447;];
tau = t9;
