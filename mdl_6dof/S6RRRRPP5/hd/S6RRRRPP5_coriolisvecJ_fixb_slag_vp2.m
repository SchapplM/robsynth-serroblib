% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 21:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 21:05:27
% EndTime: 2019-03-09 21:06:13
% DurationCPUTime: 24.95s
% Computational Cost: add. (9097->712), mult. (22550->907), div. (0->0), fcn. (14959->6), ass. (0->312)
t263 = sin(qJ(3));
t266 = cos(qJ(3));
t264 = sin(qJ(2));
t336 = qJD(1) * t264;
t316 = t266 * t336;
t210 = qJD(2) * t263 + t316;
t262 = sin(qJ(4));
t265 = cos(qJ(4));
t317 = t263 * t336;
t333 = qJD(2) * t266;
t279 = t317 - t333;
t154 = t262 * t210 + t265 * t279;
t257 = pkin(7) * t336;
t229 = -qJD(2) * pkin(2) + t257;
t179 = pkin(3) * t279 + t229;
t272 = t265 * t210 - t262 * t279;
t326 = qJD(1) * qJD(2);
t309 = t264 * t326;
t330 = qJD(3) * t264;
t267 = cos(qJ(2));
t332 = qJD(2) * t267;
t277 = -t263 * t330 + t266 * t332;
t325 = qJD(2) * qJD(3);
t173 = qJD(1) * t277 + t266 * t325;
t313 = t263 * t332;
t329 = qJD(3) * t266;
t276 = t264 * t329 + t313;
t174 = -qJD(1) * t276 - t263 * t325;
t68 = -qJD(4) * t154 + t265 * t173 + t262 * t174;
t69 = qJD(4) * t272 + t262 * t173 - t265 * t174;
t283 = Ifges(7,5) * t68 + Ifges(7,6) * t69 - Ifges(7,3) * t309;
t335 = qJD(1) * t267;
t246 = qJD(3) - t335;
t233 = -qJD(4) - t246;
t390 = pkin(4) + pkin(5);
t408 = qJ(6) * t272;
t223 = -pkin(2) * t267 - pkin(8) * t264 - pkin(1);
t202 = t223 * qJD(1);
t258 = pkin(7) * t335;
t230 = qJD(2) * pkin(8) + t258;
t160 = t266 * t202 - t230 * t263;
t126 = -pkin(9) * t210 + t160;
t114 = pkin(3) * t246 + t126;
t161 = t263 * t202 + t266 * t230;
t127 = -pkin(9) * t279 + t161;
t349 = t262 * t127;
t41 = t265 * t114 - t349;
t31 = t41 + t408;
t406 = qJD(5) - t31;
t29 = t233 * t390 + t406;
t222 = t233 * qJ(5);
t345 = t265 * t127;
t42 = t262 * t114 + t345;
t442 = qJ(6) * t154;
t32 = t42 + t442;
t30 = -t222 + t32;
t271 = -qJ(5) * t272 + t179;
t35 = -t154 * t390 + qJD(6) - t271;
t405 = qJD(5) - t41;
t37 = pkin(4) * t233 + t405;
t371 = t233 / 0.2e1;
t372 = -t233 / 0.2e1;
t38 = -t222 + t42;
t380 = t272 / 0.2e1;
t384 = -t154 / 0.2e1;
t425 = -Ifges(5,6) + Ifges(6,6);
t427 = Ifges(6,2) + Ifges(5,3);
t429 = Ifges(6,4) + Ifges(5,5);
t404 = t427 * t309 + t425 * t69 + t429 * t68;
t426 = Ifges(7,2) + Ifges(6,3);
t327 = qJD(4) * t265;
t328 = qJD(4) * t262;
t296 = pkin(2) * t264 - pkin(8) * t267;
t219 = t296 * qJD(2);
t203 = qJD(1) * t219;
t302 = pkin(7) * t309;
t102 = -qJD(3) * t161 + t266 * t203 + t263 * t302;
t57 = pkin(3) * t309 - pkin(9) * t173 + t102;
t331 = qJD(3) * t263;
t101 = t202 * t329 + t263 * t203 - t230 * t331 - t266 * t302;
t70 = pkin(9) * t174 + t101;
t11 = t114 * t327 - t127 * t328 + t262 * t57 + t265 * t70;
t12 = -t114 * t328 - t127 * t327 - t262 * t70 + t265 * t57;
t6 = qJ(5) * t309 - t233 * qJD(5) + t11;
t2 = qJ(6) * t69 + qJD(6) * t154 + t6;
t320 = t390 * t264;
t301 = qJD(2) * t320;
t3 = -qJ(6) * t68 - qJD(1) * t301 - qJD(6) * t272 - t12;
t7 = -pkin(4) * t309 - t12;
t435 = -t12 * mrSges(5,1) + t7 * mrSges(6,1) + t3 * mrSges(7,1) + t11 * mrSges(5,2) - t2 * mrSges(7,2) - t6 * mrSges(6,3);
t428 = Ifges(7,4) + Ifges(6,5);
t439 = t154 * t428;
t58 = t154 * pkin(4) + t271;
t357 = Ifges(5,4) * t272;
t81 = -Ifges(5,2) * t154 - Ifges(5,6) * t233 + t357;
t462 = -t283 + t404 - t435 + t81 * t380 + (-Ifges(7,5) * t154 + Ifges(7,6) * t272) * t372 + (-t154 * t29 - t272 * t30) * mrSges(7,3) + (t154 * t37 + t272 * t38) * mrSges(6,2) + (-t154 * t429 + t425 * t272) * t371 - t179 * (mrSges(5,1) * t272 - mrSges(5,2) * t154) - t35 * (-mrSges(7,1) * t272 - mrSges(7,2) * t154) - t58 * (mrSges(6,1) * t272 + mrSges(6,3) * t154) + (t426 * t272 - t439) * t384;
t461 = -Ifges(3,6) / 0.2e1;
t401 = -Ifges(5,4) + t428;
t460 = -Ifges(7,5) + t429;
t391 = t69 / 0.2e1;
t459 = t426 * t391;
t211 = t262 * t263 - t265 * t266;
t403 = qJD(3) + qJD(4);
t164 = t403 * t211;
t278 = t211 * t267;
t184 = qJD(1) * t278;
t455 = t164 - t184;
t393 = t68 / 0.2e1;
t381 = -t272 / 0.2e1;
t458 = t309 / 0.2e1;
t457 = qJD(2) * t461;
t456 = -qJD(1) / 0.2e1;
t424 = Ifges(6,6) - Ifges(7,6);
t212 = t262 * t266 + t263 * t265;
t165 = t403 * t212;
t183 = t212 * t335;
t339 = t165 - t183;
t366 = pkin(3) * t263;
t205 = t335 * t366 + t258;
t454 = pkin(3) * t331 - t205;
t148 = Ifges(5,4) * t154;
t453 = -Ifges(5,2) * t272 - t148;
t402 = Ifges(5,1) + Ifges(6,1) + Ifges(7,1);
t451 = t279 / 0.2e1;
t398 = -t233 * t460 + t402 * t272 - t148 + t439;
t450 = t398 / 0.2e1;
t449 = t428 * t393 + t424 * t458 + t459;
t448 = Ifges(4,6) * t279;
t285 = t160 * t266 + t161 * t263;
t447 = t285 * mrSges(4,3);
t389 = -pkin(9) - pkin(8);
t318 = qJD(3) * t389;
t217 = t263 * t318;
t218 = t266 * t318;
t231 = t389 * t263;
t232 = t389 * t266;
t110 = t265 * t217 + t262 * t218 + t231 * t327 + t232 * t328;
t216 = t296 * qJD(1);
t175 = pkin(7) * t317 + t266 * t216;
t344 = t266 * t267;
t284 = pkin(3) * t264 - pkin(9) * t344;
t140 = qJD(1) * t284 + t175;
t196 = t263 * t216;
t346 = t264 * t266;
t347 = t263 * t267;
t157 = t196 + (-pkin(7) * t346 - pkin(9) * t347) * qJD(1);
t86 = t262 * t140 + t265 * t157;
t73 = qJ(5) * t336 + t86;
t446 = t110 - t73;
t172 = t262 * t231 - t265 * t232;
t111 = qJD(4) * t172 + t217 * t262 - t265 * t218;
t85 = t140 * t265 - t262 * t157;
t445 = -t111 - t85;
t54 = t265 * t126 - t349;
t410 = pkin(3) * t327 + qJD(5) - t54;
t444 = Ifges(3,5) * qJD(2);
t350 = qJ(5) * t154;
t441 = t428 * t272;
t440 = qJ(5) * t455 - qJD(5) * t212 + t454;
t419 = t426 * t154 - t424 * t233 + t441;
t438 = -t154 * t402 + t419;
t437 = -t102 * mrSges(4,1) + t101 * mrSges(4,2);
t352 = Ifges(4,6) * t263;
t288 = Ifges(4,5) * t266 - t352;
t359 = Ifges(4,4) * t263;
t292 = Ifges(4,1) * t266 - t359;
t370 = t246 / 0.2e1;
t375 = t210 / 0.2e1;
t436 = t288 * t370 + t292 * t375;
t434 = t160 * mrSges(4,1) + t30 * mrSges(7,2) + t38 * mrSges(6,3) + t41 * mrSges(5,1) + t457 + (Ifges(3,4) * t264 + t267 * Ifges(3,2)) * t456 + Ifges(7,5) * t381 + Ifges(7,6) * t384 + Ifges(7,3) * t372 - t161 * mrSges(4,2) - t29 * mrSges(7,1) - t37 * mrSges(6,1) - t42 * mrSges(5,2);
t432 = -Ifges(5,6) / 0.2e1;
t379 = t173 / 0.2e1;
t378 = t174 / 0.2e1;
t430 = pkin(4) * t272;
t422 = -t339 * t390 - t440;
t421 = t455 * qJ(6) + qJD(1) * t320 - qJD(6) * t212 - t445;
t420 = t339 * qJ(6) + qJD(6) * t211 + t446;
t411 = -t408 + t410;
t409 = t339 * pkin(4) + t440;
t189 = t212 * t264;
t407 = t272 * t390;
t248 = pkin(7) * t344;
t182 = t263 * t223 + t248;
t399 = t309 * t460 + t401 * t69 + t402 * t68;
t360 = Ifges(4,4) * t210;
t138 = -Ifges(4,2) * t279 + Ifges(4,6) * t246 + t360;
t206 = Ifges(4,4) * t279;
t139 = t210 * Ifges(4,1) + t246 * Ifges(4,5) - t206;
t293 = mrSges(4,1) * t263 + mrSges(4,2) * t266;
t369 = t266 / 0.2e1;
t395 = -t447 + t229 * t293 - t263 * t138 / 0.2e1 + t139 * t369 + t436;
t394 = -t68 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t391 + t309 * t432;
t392 = -t69 / 0.2e1;
t386 = Ifges(4,1) * t379 + Ifges(4,4) * t378 + Ifges(4,5) * t458;
t383 = t154 / 0.2e1;
t367 = pkin(3) * t210;
t365 = pkin(7) * t263;
t364 = t68 * mrSges(7,3);
t47 = -mrSges(6,2) * t69 + mrSges(6,3) * t309;
t51 = mrSges(7,2) * t309 + mrSges(7,3) * t69;
t363 = t47 + t51;
t362 = mrSges(5,3) * t154;
t361 = mrSges(5,3) * t272;
t358 = Ifges(4,4) * t266;
t355 = Ifges(4,5) * t210;
t354 = Ifges(4,5) * t263;
t351 = Ifges(4,3) * t246;
t348 = t263 * t264;
t128 = -mrSges(6,2) * t154 - mrSges(6,3) * t233;
t129 = -mrSges(7,2) * t233 + mrSges(7,3) * t154;
t343 = t128 + t129;
t130 = mrSges(5,2) * t233 - t362;
t342 = t128 + t130;
t132 = -mrSges(5,1) * t233 - t361;
t133 = mrSges(6,1) * t233 + mrSges(6,2) * t272;
t341 = t132 - t133;
t209 = t266 * t223;
t159 = -pkin(9) * t346 + t209 + (-pkin(3) - t365) * t267;
t167 = -pkin(9) * t348 + t182;
t97 = t262 * t159 + t265 * t167;
t338 = -qJD(2) * mrSges(3,1) + mrSges(4,1) * t279 + t210 * mrSges(4,2) + mrSges(3,3) * t336;
t334 = qJD(2) * t264;
t337 = t266 * t219 + t334 * t365;
t220 = pkin(3) * t348 + t264 * pkin(7);
t321 = pkin(3) * t328;
t259 = pkin(7) * t332;
t319 = Ifges(4,5) * t173 + Ifges(4,6) * t174 + Ifges(4,3) * t309;
t180 = pkin(3) * t276 + t259;
t255 = -pkin(3) * t266 - pkin(2);
t254 = -pkin(3) * t265 - pkin(4);
t26 = -t69 * mrSges(7,1) + t68 * mrSges(7,2);
t308 = t267 * t326;
t305 = t332 / 0.2e1;
t304 = -t330 / 0.2e1;
t152 = -pkin(3) * t174 + pkin(7) * t308;
t53 = t126 * t262 + t345;
t96 = t159 * t265 - t262 * t167;
t171 = -t265 * t231 - t232 * t262;
t299 = -Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1 + Ifges(7,5) / 0.2e1;
t298 = Ifges(6,6) / 0.2e1 - Ifges(7,6) / 0.2e1 + t432;
t88 = -qJ(5) * t267 + t97;
t190 = t211 * t264;
t294 = -qJ(5) * t190 - t220;
t89 = t267 * pkin(4) - t96;
t291 = Ifges(4,1) * t263 + t358;
t290 = -Ifges(4,2) * t263 + t358;
t289 = Ifges(4,2) * t266 + t359;
t287 = -t367 - t350;
t286 = t101 * t266 - t102 * t263;
t124 = t263 * t219 + t223 * t329 + (-t264 * t333 - t267 * t331) * pkin(7);
t100 = -pkin(9) * t276 + t124;
t95 = t284 * qJD(2) + (-t248 + (pkin(9) * t264 - t223) * t263) * qJD(3) + t337;
t18 = -t262 * t100 - t159 * t328 - t167 * t327 + t265 * t95;
t282 = qJ(5) * t212 - t255;
t281 = t263 * t290;
t280 = t266 * t290;
t17 = t265 * t100 + t159 * t327 - t167 * t328 + t262 * t95;
t103 = -qJD(2) * t278 - t189 * t403;
t275 = qJ(5) * t103 - qJD(5) * t190 - t180;
t274 = qJ(5) * t68 + qJD(5) * t272 - t152;
t15 = qJ(5) * t334 - qJD(5) * t267 + t17;
t256 = Ifges(3,4) * t335;
t250 = pkin(3) * t262 + qJ(5);
t249 = -pkin(5) + t254;
t227 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t335;
t194 = Ifges(3,1) * t336 + t256 + t444;
t181 = -pkin(7) * t347 + t209;
t178 = mrSges(4,1) * t246 - mrSges(4,3) * t210;
t177 = -t246 * mrSges(4,2) - mrSges(4,3) * t279;
t176 = -pkin(7) * t316 + t196;
t151 = -mrSges(4,2) * t309 + mrSges(4,3) * t174;
t150 = mrSges(4,1) * t309 - mrSges(4,3) * t173;
t149 = pkin(4) * t211 - t282;
t137 = t351 + t355 - t448;
t136 = qJ(6) * t211 + t172;
t135 = -qJ(6) * t212 + t171;
t131 = mrSges(7,1) * t233 - mrSges(7,3) * t272;
t125 = -qJD(3) * t182 + t337;
t123 = -t211 * t390 + t282;
t116 = pkin(4) * t189 - t294;
t115 = -mrSges(4,1) * t174 + mrSges(4,2) * t173;
t105 = t173 * Ifges(4,4) + t174 * Ifges(4,2) + Ifges(4,6) * t309;
t104 = -t328 * t348 + (t346 * t403 + t313) * t265 + t277 * t262;
t94 = mrSges(5,1) * t154 + mrSges(5,2) * t272;
t93 = -mrSges(7,1) * t154 + mrSges(7,2) * t272;
t92 = mrSges(6,1) * t154 - mrSges(6,3) * t272;
t91 = t350 + t430;
t90 = -t189 * t390 + t294;
t80 = Ifges(6,4) * t272 - t233 * Ifges(6,2) + t154 * Ifges(6,6);
t78 = Ifges(5,5) * t272 - t154 * Ifges(5,6) - t233 * Ifges(5,3);
t74 = -pkin(4) * t336 - t85;
t72 = -t287 + t430;
t65 = qJ(6) * t189 + t88;
t60 = t68 * mrSges(6,2);
t56 = pkin(5) * t267 + qJ(6) * t190 + t89;
t52 = -mrSges(5,2) * t309 - mrSges(5,3) * t69;
t50 = -mrSges(6,1) * t309 + t60;
t49 = mrSges(5,1) * t309 - mrSges(5,3) * t68;
t48 = -mrSges(7,1) * t309 - t364;
t46 = -t350 - t407;
t36 = t287 - t407;
t33 = t53 + t442;
t28 = pkin(4) * t104 - t275;
t27 = mrSges(5,1) * t69 + mrSges(5,2) * t68;
t25 = mrSges(6,1) * t69 - mrSges(6,3) * t68;
t16 = -pkin(4) * t334 - t18;
t14 = -t104 * t390 + t275;
t13 = pkin(4) * t69 - t274;
t9 = qJ(6) * t104 + qJD(6) * t189 + t15;
t8 = -qJ(6) * t103 + qJD(6) * t190 - t18 - t301;
t4 = -t390 * t69 + t274;
t1 = [(-t319 / 0.2e1 - t404 / 0.2e1 - t460 * t393 - t424 * t391 - Ifges(5,6) * t392 - Ifges(4,6) * t378 - Ifges(4,5) * t379 + t283 / 0.2e1 + (t338 * pkin(7) - t279 * t290 / 0.2e1 + t444 / 0.2e1 + t436) * qJD(2) + t435 + t437) * t267 + ((-Ifges(7,5) * t190 + Ifges(7,6) * t189 + Ifges(7,3) * t267) * t456 - Ifges(7,3) * t371 + t424 * t383 + t427 * t372 + Ifges(5,6) * t384 + t434 + t460 * t380) * t334 + m(4) * (t101 * t182 + t102 * t181 + t124 * t161 + t125 * t160 + (t229 + t257) * t259) + (-t313 / 0.2e1 + t266 * t304) * t138 + t229 * t276 * mrSges(4,1) + t229 * t277 * mrSges(4,2) + (t263 * t304 + t266 * t305) * t139 + (t152 * mrSges(5,1) + t13 * mrSges(6,1) - t4 * mrSges(7,1) - t6 * mrSges(6,2) - t11 * mrSges(5,3) + t2 * mrSges(7,3) - Ifges(5,2) * t392 + t401 * t393 + t394 + t449 + t459) * t189 + (t401 * t380 - t42 * mrSges(5,3) - t38 * mrSges(6,2) + t30 * mrSges(7,3) + Ifges(7,6) * t371 + t426 * t383 + t425 * t372 + t419 / 0.2e1 - Ifges(5,2) * t384 - mrSges(7,1) * t35 + mrSges(6,1) * t58 - t81 / 0.2e1 + mrSges(5,1) * t179) * t104 + m(7) * (t14 * t35 + t2 * t65 + t29 * t8 + t3 * t56 + t30 * t9 + t4 * t90) + m(6) * (t116 * t13 + t15 * t38 + t16 * t37 + t28 * t58 + t6 * t88 + t7 * t89) + m(5) * (t11 * t97 + t12 * t96 + t152 * t220 + t17 * t42 + t179 * t180 + t18 * t41) + (-0.2e1 * pkin(1) * (mrSges(3,1) * t264 + mrSges(3,2) * t267) + (0.3e1 / 0.2e1 * t267 ^ 2 - 0.3e1 / 0.2e1 * t264 ^ 2) * Ifges(3,4)) * t326 - t105 * t348 / 0.2e1 + (t3 * mrSges(7,3) - t7 * mrSges(6,2) + t12 * mrSges(5,3) - Ifges(5,4) * t392 - t4 * mrSges(7,2) + t13 * mrSges(6,3) - t399 / 0.2e1 - t152 * mrSges(5,2) - t402 * t393 - t428 * t391) * t190 - t332 * t447 + ((0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t308 + t290 * t378 + t292 * t379 + t115 * pkin(7) + (-t101 * t263 - t102 * t266 + (t160 * t263 - t161 * t266) * qJD(3)) * mrSges(4,3) + ((t293 * t335 - t227) * pkin(7) + Ifges(4,3) * t370 + Ifges(4,5) * t375 - t448 / 0.2e1 + t457) * qJD(2)) * t264 + t194 * t305 + (t179 * mrSges(5,2) + t37 * mrSges(6,2) + t35 * mrSges(7,2) - t41 * mrSges(5,3) - t58 * mrSges(6,3) - t29 * mrSges(7,3) + Ifges(5,4) * t384 + Ifges(7,5) * t371 + t429 * t372 + t402 * t380 + t428 * t383 + t450) * t103 + ((-Ifges(4,6) * t266 - t354) * t370 - t291 * t375 + t289 * t451) * t330 + t346 * t386 + t56 * t48 + t65 * t51 + t88 * t47 + t89 * t50 + t90 * t26 + t28 * t92 + t14 * t93 + t96 * t49 + t97 * t52 + t116 * t25 + t15 * t128 + t9 * t129 + t17 * t130 + t8 * t131 + t18 * t132 + t16 * t133 + (t137 + t80 + t78 + (t264 * t288 - t429 * t190 + t425 * t189 + (-Ifges(4,3) - t427) * t267) * qJD(1)) * t334 / 0.2e1 + t124 * t177 + t125 * t178 + t180 * t94 + t181 * t150 + t182 * t151 + t220 * t27; (-t211 * t6 + t212 * t7 - t339 * t38 - t37 * t455) * mrSges(6,2) + (-mrSges(7,1) * t339 - mrSges(7,2) * t455) * t35 + (mrSges(6,1) * t339 + mrSges(6,3) * t455) * t58 + (t2 * t211 - t212 * t3 + t29 * t455 + t30 * t339) * mrSges(7,3) + (-t11 * t211 - t12 * t212 - t339 * t42 + t41 * t455) * mrSges(5,3) + (mrSges(5,1) * t339 - mrSges(5,2) * t455) * t179 + (t11 * t172 - t12 * t171 + t152 * t255 + (t110 - t86) * t42 + t445 * t41 + t454 * t179) * m(5) + t421 * t131 + (t123 * t4 + t135 * t3 + t136 * t2 + t29 * t421 + t30 * t420 + t35 * t422) * m(7) + t422 * t93 + t399 * t212 / 0.2e1 + t289 * t378 + t291 * t379 + (t94 * t366 + qJD(2) * t280 / 0.2e1 + t395) * qJD(3) + (t211 * t401 + t212 * t402) * t393 + t409 * t92 + t105 * t369 + (-t150 * t263 + t151 * t266 + (-t177 * t263 - t178 * t266) * qJD(3) + m(4) * (-qJD(3) * t285 + t286)) * pkin(8) - t341 * t111 + t342 * t110 - m(4) * (t160 * t175 + t161 * t176) + (t47 + t52) * t172 + (t50 - t49) * t171 + t286 * mrSges(4,3) + (t211 * t426 + t212 * t428) * t391 + (-t81 + t419) * (t165 / 0.2e1 - t183 / 0.2e1) + t420 * t129 + (t13 * t149 + t171 * t7 + t172 * t6 + t409 * t58 + t446 * t38 + (t111 - t74) * t37) * m(6) + t211 * t449 + t263 * t386 + (Ifges(5,4) * t212 - Ifges(5,2) * t211) * t392 + t211 * t394 + ((-t256 / 0.2e1 + pkin(1) * mrSges(3,2) * qJD(1) - t194 / 0.2e1 + (Ifges(3,5) / 0.2e1 - t280 / 0.2e1) * qJD(2) + (-m(4) * t229 + (-m(4) * pkin(2) - mrSges(4,1) * t266 + mrSges(4,2) * t263 - mrSges(3,1)) * qJD(2) - t338) * pkin(7) - t395) * t267 + (t299 * t272 + (pkin(7) * mrSges(3,2) + t354 / 0.2e1 + t461 - t299 * t212 + t298 * t211) * qJD(2) - t434 - t137 / 0.2e1 + pkin(7) * t227 + ((t352 / 0.2e1 + Ifges(3,4) / 0.2e1) * t264 + pkin(1) * mrSges(3,1)) * qJD(1) + (t281 / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t335 - t80 / 0.2e1 - t78 / 0.2e1 - t351 / 0.2e1 + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t233 - t355 / 0.2e1 - t298 * t154 - qJD(3) * t281 / 0.2e1) * t264) * qJD(1) - pkin(2) * t115 + t123 * t26 - t73 * t128 - t86 * t130 - t85 * t132 - t74 * t133 + t135 * t48 + t136 * t51 + (-t164 * t402 + t165 * t401) * t380 + (t183 * t401 - t184 * t402) * t381 + t398 * (-t164 / 0.2e1 + t184 / 0.2e1) + (-Ifges(5,4) * t164 - Ifges(5,2) * t165 + t183 * t426 - t184 * t428) * t384 + (-Ifges(5,4) * t184 - Ifges(5,2) * t183 - t164 * t428 + t165 * t426) * t383 + (-Ifges(7,5) * t164 + Ifges(7,6) * t165 + t183 * t425 - t184 * t429) * t371 + (-Ifges(7,5) * t184 + Ifges(7,6) * t183 - t164 * t429 + t165 * t425) * t372 + t149 * t25 - t176 * t177 - t175 * t178 - t205 * t94 + t4 * (-mrSges(7,1) * t211 + mrSges(7,2) * t212) + t13 * (mrSges(6,1) * t211 - mrSges(6,3) * t212) + t152 * (mrSges(5,1) * t211 + mrSges(5,2) * t212) + t255 * t27; t462 - t437 + t319 + (t272 * t401 + t438) * t381 + t363 * t250 + t410 * t128 + (t250 * t6 + t254 * t7 - t58 * t72 + t410 * t38 + (t321 - t53) * t37) * m(6) + t411 * t129 + (t2 * t250 + t249 * t3 - t35 * t36 + t411 * t30 + (t321 - t33) * t29) * m(7) + t138 * t375 + (-t210 * t94 + t262 * t52 + t265 * t49 + (t130 * t265 + (t131 - t341) * t262) * qJD(4)) * pkin(3) + t341 * t53 - t210 * (-Ifges(4,1) * t279 - t360) / 0.2e1 + t453 * t383 - t229 * (t210 * mrSges(4,1) - mrSges(4,2) * t279) - t246 * (-Ifges(4,5) * t279 - Ifges(4,6) * t210) / 0.2e1 + (-t160 * t279 + t161 * t210) * mrSges(4,3) + t154 * t450 + (-Ifges(4,2) * t210 + t139 - t206) * t451 - t72 * t92 - t36 * t93 - t54 * t130 - t33 * t131 + (-t154 * t41 + t272 * t42) * mrSges(5,3) - t160 * t177 + t161 * t178 + ((t11 * t262 + t12 * t265 + (-t262 * t41 + t265 * t42) * qJD(4)) * pkin(3) - t179 * t367 + t41 * t53 - t42 * t54) * m(5) + t249 * t48 + t254 * t50; (t341 + t361) * t42 + (-t342 - t362) * t41 + t363 * qJ(5) + (t398 + t453) * t383 + t343 * qJD(5) + (-pkin(4) * t7 + qJ(5) * t6 - t37 * t42 + t38 * t405 - t58 * t91) * m(6) + (-t357 + t438 + t441) * t381 - pkin(4) * t50 - t91 * t92 - t46 * t93 - t31 * t129 - t32 * t131 - t390 * t48 + (t2 * qJ(5) - t29 * t32 - t3 * t390 + t30 * t406 - t35 * t46) * m(7) + t462; -t364 + t60 + t343 * t233 + (t92 - t93) * t272 + (-mrSges(6,1) - mrSges(7,1)) * t309 + (t233 * t30 - t272 * t35 + t3) * m(7) + (t233 * t38 + t272 * t58 + t7) * m(6); -t154 * t129 + t272 * t131 + 0.2e1 * (t4 / 0.2e1 + t30 * t384 + t29 * t380) * m(7) + t26;];
tauc  = t1(:);
