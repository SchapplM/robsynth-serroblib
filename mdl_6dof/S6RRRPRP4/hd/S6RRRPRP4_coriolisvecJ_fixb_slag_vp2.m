% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 16:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:43:36
% EndTime: 2019-03-09 16:43:59
% DurationCPUTime: 11.54s
% Computational Cost: add. (8190->578), mult. (19974->714), div. (0->0), fcn. (13327->6), ass. (0->273)
t427 = Ifges(6,1) + Ifges(7,1);
t447 = -Ifges(6,4) + Ifges(7,5);
t426 = Ifges(7,4) + Ifges(6,5);
t446 = Ifges(6,6) - Ifges(7,6);
t424 = -Ifges(6,3) - Ifges(7,2);
t256 = sin(qJ(3));
t257 = sin(qJ(2));
t378 = cos(qJ(3));
t379 = cos(qJ(2));
t223 = t256 * t379 + t257 * t378;
t254 = qJD(2) + qJD(3);
t184 = t254 * t223;
t170 = t184 * qJD(1);
t255 = sin(qJ(5));
t298 = t378 * t379;
t278 = qJD(1) * t298;
t331 = qJD(1) * t257;
t211 = t256 * t331 - t278;
t258 = cos(qJ(5));
t279 = t258 * t211 - t254 * t255;
t435 = qJD(5) * t279;
t99 = t170 * t255 + t435;
t399 = t99 / 0.2e1;
t188 = t211 * t255 + t254 * t258;
t100 = qJD(5) * t188 - t258 * t170;
t397 = -t100 / 0.2e1;
t335 = t256 * t257;
t293 = t254 * t335;
t169 = qJD(1) * t293 - t254 * t278;
t395 = -t169 / 0.2e1;
t445 = mrSges(5,2) - mrSges(4,1);
t398 = pkin(3) + pkin(9);
t232 = (-pkin(8) - pkin(7)) * t257;
t225 = qJD(1) * t232;
t219 = qJD(2) * pkin(2) + t225;
t200 = t378 * t219;
t321 = t379 * pkin(7);
t233 = pkin(8) * t379 + t321;
t226 = t233 * qJD(1);
t336 = t256 * t226;
t177 = -t200 + t336;
t212 = t223 * qJD(1);
t373 = t212 * pkin(4);
t276 = t177 + t373;
t434 = qJD(4) + t276;
t102 = -t254 * t398 + t434;
t246 = -pkin(2) * t379 - pkin(1);
t231 = qJD(1) * t246;
t265 = -t212 * qJ(4) + t231;
t104 = t211 * t398 + t265;
t28 = t102 * t258 - t104 * t255;
t29 = t102 * t255 + t104 * t258;
t282 = t255 * t28 - t258 * t29;
t325 = qJD(1) * qJD(2);
t306 = t257 * t325;
t240 = pkin(2) * t306;
t272 = qJ(4) * t169 - qJD(4) * t212 + t240;
t26 = t170 * t398 + t272;
t348 = qJD(5) * t29;
t315 = t378 * t226;
t178 = t256 * t219 + t315;
t266 = qJD(2) * t233;
t217 = t378 * t266;
t227 = qJD(2) * t232;
t220 = qJD(1) * t227;
t92 = qJD(1) * t217 + qJD(3) * t178 + t256 * t220;
t45 = -t169 * pkin(4) + t92;
t5 = -t255 * t26 + t258 * t45 - t348;
t3 = pkin(5) * t169 - t5;
t372 = t258 * t3;
t208 = qJD(5) + t212;
t416 = qJD(6) - t28;
t21 = -pkin(5) * t208 + t416;
t22 = qJ(6) * t208 + t29;
t284 = t21 * t255 + t22 * t258;
t326 = qJD(5) * t258;
t327 = qJD(5) * t255;
t4 = t102 * t326 - t104 * t327 + t255 * t45 + t258 * t26;
t1 = -qJ(6) * t169 + qJD(6) * t208 + t4;
t374 = t1 * t255;
t402 = -qJD(5) * t284 - t374;
t371 = t4 * t255;
t406 = t5 * t258 + t371;
t261 = m(6) * (-qJD(5) * t282 + t406) + m(7) * (-t372 - t402);
t49 = -mrSges(6,1) * t169 - mrSges(6,3) * t99;
t50 = t169 * mrSges(7,1) + t99 * mrSges(7,2);
t367 = t49 - t50;
t48 = -mrSges(7,2) * t100 - mrSges(7,3) * t169;
t51 = mrSges(6,2) * t169 - mrSges(6,3) * t100;
t368 = t48 + t51;
t444 = t368 * t255 + t367 * t258 + t261;
t443 = -Ifges(4,5) + Ifges(5,4);
t442 = -Ifges(4,6) + Ifges(5,5);
t441 = t447 * t100 - t426 * t169 + t427 * t99;
t186 = Ifges(6,4) * t279;
t360 = Ifges(7,5) * t279;
t423 = t188 * t427 + t426 * t208 + t186 - t360;
t316 = t378 * t225;
t180 = t316 - t336;
t307 = qJD(3) * t378;
t299 = pkin(2) * t307;
t236 = t299 + qJD(4);
t440 = t180 - t236;
t277 = -pkin(5) * t258 - qJ(6) * t255 - pkin(4);
t439 = pkin(5) * t326 + qJ(6) * t327 - t258 * qJD(6) - t212 * t277 + qJD(4) + t336;
t438 = t255 * t426 + t258 * t446;
t358 = Ifges(7,5) * t258;
t361 = Ifges(6,4) * t258;
t437 = t255 * t427 - t358 + t361;
t179 = t225 * t256 + t315;
t329 = qJD(3) * t256;
t319 = pkin(2) * t329;
t436 = t319 - t179;
t206 = Ifges(4,4) * t211;
t433 = t212 * Ifges(4,1) + t254 * Ifges(4,5) + t188 * t426 - t208 * t424 + t279 * t446 - t206;
t432 = -t99 * Ifges(7,5) / 0.2e1 + t169 * Ifges(7,6) / 0.2e1 + Ifges(6,4) * t399 + Ifges(6,6) * t395 + (Ifges(7,3) + Ifges(6,2)) * t397;
t376 = pkin(4) * t211;
t128 = t178 - t376;
t252 = t254 * qJ(4);
t116 = t128 + t252;
t291 = mrSges(7,1) * t258 + mrSges(7,3) * t255;
t292 = mrSges(6,1) * t258 - mrSges(6,2) * t255;
t42 = -pkin(5) * t279 - qJ(6) * t188 + t116;
t363 = Ifges(6,4) * t188;
t85 = Ifges(6,2) * t279 + Ifges(6,6) * t208 + t363;
t431 = t116 * t292 + t291 * t42 - t258 * t85 / 0.2e1;
t139 = t211 * pkin(3) + t265;
t430 = t231 * mrSges(4,1) - t139 * mrSges(5,2);
t405 = -qJD(4) - t177;
t141 = -pkin(3) * t254 - t405;
t428 = t141 * mrSges(5,1) + t28 * mrSges(6,1) - t21 * mrSges(7,1) + t231 * mrSges(4,2) - t29 * mrSges(6,2) + t177 * mrSges(4,3) - t139 * mrSges(5,3) + t22 * mrSges(7,3);
t421 = -t200 + t439;
t420 = -t316 + t299 + t439;
t222 = -t298 + t335;
t270 = -t223 * qJ(4) + t246;
t129 = t222 * t398 + t270;
t189 = -t378 * t232 + t233 * t256;
t155 = pkin(4) * t223 + t189;
t415 = t258 * t129 + t255 * t155;
t414 = t373 - t440;
t192 = -mrSges(4,2) * t254 - mrSges(4,3) * t211;
t194 = mrSges(5,1) * t211 - mrSges(5,3) * t254;
t413 = t192 - t194;
t353 = t212 * mrSges(4,3);
t412 = t212 * mrSges(5,1) + t254 * t445 + t353;
t320 = t378 * pkin(2);
t245 = -t320 - pkin(3);
t241 = -pkin(9) + t245;
t411 = t241 * t326 + t255 * t319;
t410 = t241 * t327 - t258 * t319;
t409 = -t100 * t446 + t169 * t424 + t426 * t99;
t309 = qJD(1) * t379;
t408 = t257 * pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t309) + (qJD(2) * mrSges(3,1) - mrSges(3,3) * t331) * t321;
t407 = t255 * pkin(5) - qJ(6) * t258;
t404 = Ifges(5,6) + Ifges(4,4) / 0.2e1;
t190 = t256 * t232 + t233 * t378;
t403 = -t169 * t189 - t170 * t190 + t223 * t92;
t183 = -t254 * t298 + t293;
t330 = qJD(2) * t257;
t250 = pkin(2) * t330;
t271 = qJ(4) * t183 - qJD(4) * t223 + t250;
t37 = t184 * t398 + t271;
t106 = qJD(3) * t190 + t256 * t227 + t217;
t64 = -t183 * pkin(4) + t106;
t9 = -qJD(5) * t415 - t255 * t37 + t258 * t64;
t401 = t5 * mrSges(6,1) - t3 * mrSges(7,1) - t4 * mrSges(6,2) + t1 * mrSges(7,3);
t396 = t100 / 0.2e1;
t393 = t279 / 0.2e1;
t392 = -t279 / 0.2e1;
t391 = -t188 / 0.2e1;
t390 = t188 / 0.2e1;
t389 = -t208 / 0.2e1;
t388 = t208 / 0.2e1;
t387 = -t211 / 0.2e1;
t386 = t211 / 0.2e1;
t385 = -t212 / 0.2e1;
t384 = t212 / 0.2e1;
t382 = -t254 / 0.2e1;
t381 = t254 / 0.2e1;
t377 = pkin(2) * t256;
t375 = pkin(5) * t211;
t366 = mrSges(6,3) * t279;
t365 = mrSges(6,3) * t188;
t364 = Ifges(3,4) * t257;
t362 = Ifges(6,4) * t255;
t359 = Ifges(7,5) * t255;
t356 = t169 * mrSges(5,1);
t354 = t189 * t92;
t352 = t212 * Ifges(4,4);
t351 = t212 * Ifges(5,6);
t157 = -t252 - t178;
t347 = t157 * t212;
t344 = t184 * t255;
t343 = t184 * t258;
t342 = t212 * t255;
t341 = t212 * t258;
t340 = t222 * t255;
t172 = pkin(3) * t212 + qJ(4) * t211;
t248 = pkin(2) * t331;
t142 = t172 + t248;
t207 = t212 * pkin(9);
t112 = t142 + t207;
t130 = t179 - t376;
t41 = t258 * t112 + t255 * t130;
t115 = -mrSges(6,1) * t279 + mrSges(6,2) * t188;
t334 = t115 - t194;
t123 = t172 + t207;
t44 = t258 * t123 + t255 * t128;
t132 = mrSges(7,2) * t279 + mrSges(7,3) * t208;
t133 = -mrSges(6,2) * t208 + t366;
t333 = t132 + t133;
t134 = mrSges(6,1) * t208 - t365;
t135 = -mrSges(7,1) * t208 + mrSges(7,2) * t188;
t332 = -t134 + t135;
t318 = Ifges(3,4) * t379;
t308 = qJD(2) * t379;
t304 = -t327 / 0.2e1;
t303 = t326 / 0.2e1;
t296 = qJD(1) * t308;
t228 = qJ(4) + t407;
t288 = Ifges(6,2) * t258 + t362;
t285 = -Ifges(7,3) * t258 + t359;
t40 = -t112 * t255 + t130 * t258;
t43 = -t123 * t255 + t128 * t258;
t53 = -t129 * t255 + t155 * t258;
t275 = Ifges(3,5) * t379 - Ifges(3,6) * t257;
t274 = t222 * t326 + t344;
t273 = t222 * t327 - t343;
t8 = -t129 * t327 + t155 * t326 + t255 * t64 + t258 * t37;
t213 = t256 * t266;
t91 = -qJD(1) * t213 + t219 * t307 + t378 * t220 - t226 * t329;
t105 = -t378 * t227 - t232 * t307 + t233 * t329 + t213;
t269 = pkin(1) * (mrSges(3,1) * t257 + mrSges(3,2) * t379);
t268 = t257 * (Ifges(3,1) * t379 - t364);
t267 = (Ifges(3,2) * t379 + t364) * qJD(1);
t79 = -t254 * qJD(4) - t91;
t38 = -pkin(4) * t170 - t79;
t11 = pkin(5) * t100 - qJ(6) * t99 - qJD(6) * t188 + t38;
t149 = t254 * Ifges(5,5) + t211 * Ifges(5,3) - t351;
t205 = Ifges(5,6) * t211;
t150 = t254 * Ifges(5,4) - t212 * Ifges(5,2) + t205;
t151 = -t211 * Ifges(4,2) + t254 * Ifges(4,6) + t352;
t185 = Ifges(7,5) * t188;
t82 = Ifges(7,6) * t208 - Ifges(7,3) * t279 + t185;
t260 = (t303 + t341 / 0.2e1) * t82 + t431 * qJD(5) + (-t352 + t149) * t385 + t358 * t396 + t361 * t397 - t91 * mrSges(4,2) - t79 * mrSges(5,3) + (-t21 * t342 - t22 * t341 + t372) * mrSges(7,2) + (t351 + t151) * t384 + (-t29 * t341 + (t327 + t342) * t28) * mrSges(6,3) + (t359 - t362) * t399 + t445 * t92 + (-Ifges(4,2) * t212 - t206 + t433) * t386 + (t205 + t150) * t387 + (t38 * mrSges(6,1) + t11 * mrSges(7,1) - Ifges(6,2) * t397 + Ifges(7,3) * t396 - t395 * t446 - t432) * t255 + (Ifges(5,3) * t387 + t285 * t393 + t288 * t392 + t382 * t442 + t389 * t438 + t391 * t437 - t430 + t431) * t212 + t442 * t170 + t443 * t169 + (-Ifges(4,1) * t385 + Ifges(5,2) * t384 - Ifges(6,6) * t392 - Ifges(7,6) * t393 + t382 * t443 + t389 * t424 - t391 * t426 + t428) * t211 + (-t288 / 0.2e1 + t285 / 0.2e1) * t435 + (-t342 / 0.2e1 + t304) * t423 - (t188 * t437 + t208 * t438) * qJD(5) / 0.2e1 + (t441 / 0.2e1 + t38 * mrSges(6,2) - t11 * mrSges(7,3) + t426 * t395 + t427 * t399) * t258;
t247 = Ifges(3,4) * t309;
t242 = qJ(4) + t377;
t218 = t228 + t377;
t210 = Ifges(3,1) * t331 + Ifges(3,5) * qJD(2) + t247;
t209 = Ifges(3,6) * qJD(2) + t267;
t204 = t211 * qJ(6);
t176 = t222 * pkin(3) + t270;
t175 = -mrSges(5,2) * t211 - mrSges(5,3) * t212;
t174 = mrSges(4,1) * t211 + mrSges(4,2) * t212;
t156 = -t222 * pkin(4) + t190;
t114 = -mrSges(7,1) * t279 - mrSges(7,3) * t188;
t113 = pkin(5) * t188 - qJ(6) * t279;
t90 = t222 * t277 + t190;
t70 = pkin(3) * t184 + t271;
t63 = -pkin(4) * t184 - t105;
t52 = pkin(3) * t170 + t272;
t47 = -pkin(5) * t223 - t53;
t46 = qJ(6) * t223 + t415;
t33 = -t43 + t375;
t32 = -t204 + t44;
t31 = -t40 + t375;
t30 = -t204 + t41;
t25 = mrSges(6,1) * t100 + mrSges(6,2) * t99;
t24 = mrSges(7,1) * t100 - mrSges(7,3) * t99;
t12 = (qJD(5) * t407 - qJD(6) * t255) * t222 + t277 * t184 - t105;
t7 = pkin(5) * t183 - t9;
t6 = -qJ(6) * t183 + qJD(6) * t223 + t8;
t2 = [(t21 * t274 - t22 * t273 + t3 * t340) * mrSges(7,2) + t174 * t250 + (t268 - 0.2e1 * t269) * t325 + (Ifges(7,5) * t274 + Ifges(7,3) * t273) * t392 + (-t273 * t446 + t274 * t426) * t388 + ((t1 * mrSges(7,2) + t4 * mrSges(6,3) + t432) * t258 + t423 * t303 + mrSges(4,1) * t240 + t79 * mrSges(5,1) - t52 * mrSges(5,2) - t91 * mrSges(4,3) - t11 * t291 + t285 * t396 + t288 * t397 - t38 * t292 + t85 * t304 + t437 * t399 + (Ifges(5,3) + Ifges(4,2)) * t170 + t404 * t169 + (-Ifges(4,4) + t438) * t395) * t222 + (qJD(1) * (Ifges(3,1) * t257 + t318) + t210) * t308 / 0.2e1 - (t267 + t209) * t330 / 0.2e1 + (-t433 / 0.2e1 + t150 / 0.2e1 + t443 * t381 - Ifges(4,1) * t384 + Ifges(5,2) * t385 + Ifges(5,6) * t386 - Ifges(4,4) * t387 - t428 - t426 * t390 + t424 * t388 - Ifges(7,6) * t392 - Ifges(6,6) * t393) * t183 + (Ifges(6,4) * t274 - Ifges(6,2) * t273) * t393 + t156 * t25 + t8 * t133 + t9 * t134 + t7 * t135 + m(7) * (t1 * t46 + t11 * t90 + t12 * t42 + t21 * t7 + t22 * t6 + t3 * t47) + t6 * t132 + t12 * t114 + t63 * t115 + t90 * t24 + t53 * t49 + t46 * t48 + t47 * t50 + m(5) * (t139 * t70 + t176 * t52 - t190 * t79 + t354) + m(4) * (t190 * t91 + 0.2e1 * t231 * t250 + t354) + t423 * t344 / 0.2e1 + (mrSges(4,2) * t240 - t52 * mrSges(5,3) - Ifges(5,2) * t169 + Ifges(6,6) * t397 + Ifges(7,6) * t396 + t426 * t399 - t404 * t170 + (Ifges(4,1) - t424) * t395 + t401) * t223 + m(6) * (t116 * t63 + t156 * t38 + t28 * t9 + t29 * t8 + t4 * t415 + t5 * t53) + t415 * t51 + (t275 * qJD(2) / 0.2e1 - t408) * qJD(2) + (-Ifges(4,1) * t169 - Ifges(4,4) * t170 + t409) * t223 / 0.2e1 + (m(4) * t177 + m(5) * t141 + t412) * t106 + (-m(4) * t178 + m(5) * t157 - t413) * t105 + (t85 / 0.2e1 - t82 / 0.2e1) * t343 + (-t273 * t29 - t274 * t28 - t340 * t5) * mrSges(6,3) + (-t178 * t184 + t403) * mrSges(4,3) + (t157 * t184 + t403) * mrSges(5,1) + (-Ifges(3,2) * t257 + t318) * t296 + t116 * (mrSges(6,1) * t273 + mrSges(6,2) * t274) + t42 * (mrSges(7,1) * t273 - mrSges(7,3) * t274) + t70 * t175 + t176 * (-mrSges(5,2) * t170 + mrSges(5,3) * t169) + (qJD(5) * t82 + t441) * t340 / 0.2e1 + (t149 / 0.2e1 - t151 / 0.2e1 - Ifges(4,4) * t384 + Ifges(5,6) * t385 + Ifges(5,3) * t386 - Ifges(4,2) * t387 + t442 * t381 + t430) * t184 + t246 * (mrSges(4,1) * t170 - mrSges(4,2) * t169) + (t273 * t447 + t427 * t274) * t390; t420 * t114 + t436 * t412 + (-t139 * t142 + t141 * t436 + t157 * t440 - t242 * t79 + t245 * t92) * m(5) - (-Ifges(3,2) * t331 + t210 + t247) * t309 / 0.2e1 + Ifges(3,5) * t296 + t192 * t299 + (-t21 * t31 - t22 * t30 + t11 * t218 + (-t21 * t258 + t22 * t255) * t319 + t420 * t42) * m(7) + (-t28 * t40 - t29 * t41 + t242 * t38 + (t255 * t29 + t258 * t28) * t319 + t414 * t116) * m(6) + t414 * t115 + (-t40 - t410) * t134 + (-t31 + t410) * t135 + (-t41 + t411) * t133 + (-t30 + t411) * t132 - t413 * t180 + t260 + t444 * t241 + (-mrSges(3,1) * t296 + mrSges(3,2) * t306) * pkin(7) + (-t29 * t326 - t406) * mrSges(6,3) - t245 * t356 + (t169 * t320 - t170 * t377) * mrSges(4,3) + (-t21 * t327 - t22 * t326 - t374) * mrSges(7,2) + t209 * t331 / 0.2e1 - t275 * t325 / 0.2e1 + (-t170 * t242 - t347) * mrSges(5,1) - t174 * t248 + (t408 + (-t268 / 0.2e1 + t269) * qJD(1)) * qJD(1) - Ifges(3,6) * t306 + (-t177 * t179 - t178 * t180 - t231 * t248 + (-t378 * t92 + t256 * t91 + (t177 * t256 + t178 * t378) * qJD(3)) * pkin(2)) * m(4) + t178 * t353 - t142 * t175 + t218 * t24 - t236 * t194 + t242 * t25; (pkin(3) * t169 - qJ(4) * t170 - t347) * mrSges(5,1) - t44 * t133 - t43 * t134 - t33 * t135 + t276 * t115 - t32 * t132 + qJ(4) * t25 + t402 * mrSges(7,2) + t260 + (-t371 + (-t5 - t348) * t258) * mrSges(6,3) + t421 * t114 - ((t255 * t332 + t258 * t333) * qJD(5) + t444) * t398 + (t353 - t412) * t178 + t413 * t177 + t334 * qJD(4) - t172 * t175 + t228 * t24 + (t11 * t228 - t21 * t33 - t22 * t32 + t42 * t421) * m(7) + (t38 * qJ(4) + t116 * t434 - t28 * t43 - t29 * t44) * m(6) + (-pkin(3) * t92 - qJ(4) * t79 - t139 * t172 - t141 * t178 + t157 * t405) * m(5); -t356 + t212 * t175 + (-t114 - t334) * t254 + (t208 * t333 + t367) * t258 + (t208 * t332 + t368) * t255 - m(6) * (t116 * t254 + t212 * t282) - m(7) * (-t212 * t284 + t254 * t42) + t261 + (t139 * t212 + t157 * t254 + t92) * m(5); qJD(6) * t132 - t113 * t114 + qJ(6) * t48 - pkin(5) * t50 + t401 + (-Ifges(6,2) * t188 + t186 + t423) * t392 + (-t332 + t365) * t29 + (-pkin(5) * t3 + qJ(6) * t1 - t113 * t42 - t21 * t29 + t22 * t416) * m(7) + (t188 * t22 - t21 * t279) * mrSges(7,2) - t42 * (mrSges(7,1) * t188 - mrSges(7,3) * t279) - t116 * (mrSges(6,1) * t188 + mrSges(6,2) * t279) + (-t188 * t446 + t279 * t426) * t389 + (t279 * t427 + t185 - t363 + t82) * t391 + (Ifges(7,3) * t188 + t360) * t393 + (-t333 + t366) * t28 + t85 * t390 + t409; t188 * t114 - t208 * t132 + 0.2e1 * (t3 / 0.2e1 + t42 * t390 + t22 * t389) * m(7) + t50;];
tauc  = t2(:);
