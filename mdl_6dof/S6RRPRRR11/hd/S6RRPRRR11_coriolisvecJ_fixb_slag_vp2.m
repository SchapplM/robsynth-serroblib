% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,d6]';
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
% Datum: 2019-03-09 14:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR11_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR11_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 14:29:17
% EndTime: 2019-03-09 14:29:44
% DurationCPUTime: 15.16s
% Computational Cost: add. (14525->688), mult. (32473->955), div. (0->0), fcn. (21062->8), ass. (0->324)
t291 = sin(qJ(4));
t295 = cos(qJ(4));
t296 = cos(qJ(2));
t359 = qJD(1) * t296;
t235 = -qJD(2) * t291 - t295 * t359;
t357 = qJD(2) * t295;
t236 = -t291 * t359 + t357;
t290 = sin(qJ(5));
t294 = cos(qJ(5));
t164 = t235 * t290 + t236 * t294;
t289 = sin(qJ(6));
t293 = cos(qJ(6));
t330 = t294 * t235 - t236 * t290;
t105 = t164 * t293 + t289 * t330;
t279 = pkin(7) * t359;
t244 = pkin(3) * t359 + t279;
t288 = qJD(2) * qJ(3);
t220 = t288 + t244;
t175 = -pkin(4) * t235 + t220;
t113 = -pkin(5) * t330 + t175;
t468 = pkin(10) * t330;
t297 = -pkin(2) - pkin(8);
t292 = sin(qJ(2));
t335 = -qJ(3) * t292 - pkin(1);
t227 = t296 * t297 + t335;
t193 = t227 * qJD(1);
t283 = t292 * qJD(1);
t277 = pkin(7) * t283;
t243 = -pkin(3) * t283 - t277;
t442 = qJD(3) - t243;
t198 = qJD(2) * t297 + t442;
t130 = -t193 * t291 + t295 * t198;
t114 = -pkin(9) * t236 + t130;
t270 = t283 + qJD(4);
t106 = pkin(4) * t270 + t114;
t131 = t193 * t295 + t198 * t291;
t115 = pkin(9) * t235 + t131;
t112 = t294 * t115;
t51 = t106 * t290 + t112;
t41 = t51 + t468;
t377 = t289 * t41;
t345 = qJD(4) + qJD(5);
t261 = t283 + t345;
t480 = pkin(10) * t164;
t110 = t290 * t115;
t50 = t294 * t106 - t110;
t40 = t50 - t480;
t39 = pkin(5) * t261 + t40;
t14 = t293 * t39 - t377;
t375 = t293 * t41;
t15 = t289 * t39 + t375;
t347 = qJD(1) * qJD(2);
t336 = t296 * t347;
t266 = Ifges(7,3) * t336;
t392 = Ifges(7,4) * t105;
t253 = qJD(6) + t261;
t407 = -t253 / 0.2e1;
t423 = t105 / 0.2e1;
t424 = -t105 / 0.2e1;
t473 = -t164 * t289 + t293 * t330;
t426 = -t473 / 0.2e1;
t350 = qJD(5) * t294;
t351 = qJD(5) * t290;
t346 = qJD(2) * qJD(4);
t352 = qJD(4) * t296;
t358 = qJD(2) * t292;
t183 = -t291 * t346 + (t291 * t358 - t295 * t352) * qJD(1);
t337 = t292 * t347;
t269 = pkin(2) * t337;
t319 = pkin(8) * t292 - qJ(3) * t296;
t355 = qJD(3) * t292;
t302 = qJD(2) * t319 - t355;
t179 = qJD(1) * t302 + t269;
t356 = qJD(2) * t296;
t430 = pkin(3) + pkin(7);
t246 = t430 * t356;
t226 = qJD(1) * t246;
t75 = -qJD(4) * t131 - t179 * t291 + t295 * t226;
t55 = pkin(4) * t336 - pkin(9) * t183 + t75;
t338 = t291 * t352;
t303 = t292 * t357 + t338;
t184 = qJD(1) * t303 - t295 * t346;
t353 = qJD(4) * t295;
t354 = qJD(4) * t291;
t74 = t295 * t179 - t193 * t354 + t198 * t353 + t291 * t226;
t61 = pkin(9) * t184 + t74;
t12 = t106 * t350 - t115 * t351 + t290 * t55 + t294 * t61;
t86 = -qJD(5) * t164 - t183 * t290 + t184 * t294;
t10 = pkin(10) * t86 + t12;
t13 = -qJD(5) * t51 - t290 * t61 + t294 * t55;
t85 = qJD(5) * t330 + t183 * t294 + t184 * t290;
t9 = pkin(5) * t336 - pkin(10) * t85 + t13;
t2 = qJD(6) * t14 + t10 * t293 + t289 * t9;
t3 = -qJD(6) * t15 - t10 * t289 + t293 * t9;
t30 = qJD(6) * t473 + t289 * t86 + t293 * t85;
t31 = -qJD(6) * t105 - t289 * t85 + t293 * t86;
t441 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t30 + Ifges(7,6) * t31;
t47 = Ifges(7,2) * t473 + Ifges(7,6) * t253 + t392;
t98 = Ifges(7,4) * t473;
t48 = Ifges(7,1) * t105 + Ifges(7,5) * t253 + t98;
t501 = (Ifges(7,1) * t473 - t392) * t424 + (Ifges(7,5) * t473 - Ifges(7,6) * t105) * t407 + (t105 * t15 + t14 * t473) * mrSges(7,3) - t113 * (mrSges(7,1) * t105 + mrSges(7,2) * t473) + t47 * t423 + t266 + t441 + (-Ifges(7,2) * t105 + t48 + t98) * t426;
t278 = pkin(2) * t283;
t202 = qJD(1) * t319 + t278;
t151 = -t202 * t291 + t295 * t244;
t309 = -pkin(9) * t291 * t292 + pkin(4) * t296;
t397 = pkin(9) - t297;
t500 = -qJD(1) * t309 + t397 * t354 - t151;
t152 = t295 * t202 + t291 * t244;
t250 = t397 * t295;
t339 = t295 * t283;
t499 = pkin(9) * t339 + qJD(4) * t250 + t152;
t157 = Ifges(6,4) * t330;
t267 = Ifges(6,3) * t336;
t393 = Ifges(6,4) * t164;
t405 = -t261 / 0.2e1;
t417 = -t164 / 0.2e1;
t419 = -t330 / 0.2e1;
t440 = t13 * mrSges(6,1) - t12 * mrSges(6,2) + Ifges(6,5) * t85 + Ifges(6,6) * t86;
t94 = Ifges(6,1) * t164 + Ifges(6,5) * t261 + t157;
t498 = t267 + t440 + (Ifges(6,5) * t330 - Ifges(6,6) * t164) * t405 + (t164 * t51 + t330 * t50) * mrSges(6,3) + (-Ifges(6,2) * t164 + t157 + t94) * t419 - t175 * (mrSges(6,1) * t164 + mrSges(6,2) * t330) + (Ifges(6,1) * t330 - t393) * t417 + t501;
t364 = t294 * t295;
t172 = -t290 * t354 - t291 * t351 + t345 * t364;
t366 = t290 * t291;
t194 = -t283 * t366 + t294 * t339;
t362 = t172 + t194;
t311 = t290 * t295 + t294 * t291;
t173 = t345 * t311;
t304 = t311 * t292;
t195 = qJD(1) * t304;
t496 = t173 + t195;
t249 = t397 * t291;
t177 = -t294 * t249 - t290 * t250;
t450 = -qJD(5) * t177 + t499 * t290 + t294 * t500;
t449 = t249 * t351 - t250 * t350 + t290 * t500 - t499 * t294;
t310 = -t364 + t366;
t166 = t289 * t310 - t293 * t311;
t301 = qJD(6) * t166 - t172 * t289 - t293 * t173;
t332 = t194 * t289 + t293 * t195;
t372 = t332 - t301;
t445 = -t289 * t311 - t293 * t310;
t125 = t194 * t293 - t195 * t289;
t66 = qJD(6) * t445 + t172 * t293 - t173 * t289;
t485 = t66 + t125;
t493 = t14 * t372 - t15 * t485 + t166 * t2 - t445 * t3;
t490 = -pkin(5) * t359 + t496 * pkin(10) + t450;
t489 = t362 * pkin(10) - t449;
t396 = Ifges(5,4) * t236;
t148 = t235 * Ifges(5,2) + t270 * Ifges(5,6) + t396;
t228 = Ifges(5,4) * t235;
t149 = t236 * Ifges(5,1) + t270 * Ifges(5,5) + t228;
t316 = t130 * t291 - t131 * t295;
t395 = Ifges(5,4) * t291;
t321 = Ifges(5,2) * t295 + t395;
t394 = Ifges(5,4) * t295;
t323 = Ifges(5,1) * t291 + t394;
t326 = mrSges(5,1) * t295 - mrSges(5,2) * t291;
t387 = Ifges(5,6) * t295;
t391 = Ifges(5,5) * t291;
t401 = -t295 / 0.2e1;
t402 = -t291 / 0.2e1;
t403 = -t270 / 0.2e1;
t410 = -t236 / 0.2e1;
t411 = -t235 / 0.2e1;
t484 = t316 * mrSges(5,3) + t148 * t401 + t149 * t402 + t220 * t326 + (t387 + t391) * t403 + t321 * t411 + t323 * t410;
t483 = qJD(3) + t277;
t479 = -mrSges(3,1) + mrSges(4,2);
t275 = pkin(4) * t294 + pkin(5);
t348 = qJD(6) * t293;
t349 = qJD(6) * t289;
t365 = t290 * t293;
t57 = -t114 * t290 - t112;
t44 = t57 - t468;
t58 = t294 * t114 - t110;
t45 = t58 - t480;
t478 = t289 * t45 - t293 * t44 - t275 * t349 + (-t290 * t348 + (-t289 * t294 - t365) * qJD(5)) * pkin(4);
t367 = t289 * t290;
t477 = -t289 * t44 - t293 * t45 + t275 * t348 + (-t290 * t349 + (t293 * t294 - t367) * qJD(5)) * pkin(4);
t340 = -pkin(4) * t295 - pkin(3);
t444 = pkin(4) * t353 - t283 * t340 + t483;
t251 = -pkin(2) * t296 + t335;
t221 = t251 * qJD(1);
t255 = -t279 - t288;
t464 = qJD(2) / 0.2e1;
t465 = -qJD(2) / 0.2e1;
t466 = -qJD(1) / 0.2e1;
t476 = Ifges(3,6) * t464 + (Ifges(3,4) * t292 + t296 * Ifges(3,2)) * qJD(1) / 0.2e1 + Ifges(4,5) * t465 + (-Ifges(4,6) * t292 - t296 * Ifges(4,3)) * t466 + t221 * mrSges(4,2) - t255 * mrSges(4,1) + t484;
t475 = -t12 * t311 + t13 * t310 - t362 * t51 + t496 * t50;
t93 = Ifges(6,2) * t330 + Ifges(6,6) * t261 + t393;
t469 = t93 / 0.2e1;
t176 = t249 * t290 - t294 * t250;
t141 = pkin(10) * t310 + t176;
t142 = -pkin(10) * t311 + t177;
t80 = t141 * t293 - t142 * t289;
t463 = qJD(6) * t80 + t289 * t490 - t489 * t293;
t81 = t141 * t289 + t142 * t293;
t462 = -qJD(6) * t81 + t489 * t289 + t293 * t490;
t446 = pkin(5) * t362 + t444;
t259 = t430 * t292;
t240 = t295 * t259;
t334 = pkin(9) * t296 - t227;
t144 = pkin(4) * t292 + t291 * t334 + t240;
t239 = t291 * t259;
t171 = t295 * t227 + t239;
t363 = t295 * t296;
t150 = -pkin(9) * t363 + t171;
t89 = t290 * t144 + t294 * t150;
t443 = t291 * t74 + t295 * t75;
t439 = t296 * t345;
t438 = t75 * mrSges(5,1) - t74 * mrSges(5,2) + Ifges(5,5) * t183 + Ifges(5,6) * t184;
t436 = t30 / 0.2e1;
t435 = t31 / 0.2e1;
t432 = t85 / 0.2e1;
t431 = t86 / 0.2e1;
t429 = pkin(1) * mrSges(3,1);
t428 = pkin(1) * mrSges(3,2);
t425 = t473 / 0.2e1;
t203 = t310 * t296;
t204 = t311 * t296;
t139 = t203 * t293 + t204 * t289;
t422 = t139 / 0.2e1;
t140 = t203 * t289 - t204 * t293;
t421 = t140 / 0.2e1;
t420 = t148 / 0.2e1;
t418 = t330 / 0.2e1;
t416 = t164 / 0.2e1;
t415 = t166 / 0.2e1;
t414 = t445 / 0.2e1;
t413 = t203 / 0.2e1;
t412 = -t204 / 0.2e1;
t409 = -t311 / 0.2e1;
t408 = -t310 / 0.2e1;
t406 = t253 / 0.2e1;
t404 = t261 / 0.2e1;
t390 = Ifges(5,5) * t295;
t389 = Ifges(4,6) * t296;
t388 = Ifges(5,6) * t291;
t371 = qJD(2) * mrSges(3,2);
t272 = t291 * pkin(4) + qJ(3);
t174 = -mrSges(5,1) * t235 + mrSges(5,2) * t236;
t257 = -mrSges(4,1) * t359 - qJD(2) * mrSges(4,3);
t360 = -t257 + t174;
t260 = t430 * t296;
t344 = Ifges(3,5) / 0.2e1 - Ifges(4,4) / 0.2e1;
t343 = -Ifges(3,6) / 0.2e1 + Ifges(4,5) / 0.2e1;
t342 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t341 = m(4) * pkin(7) + mrSges(4,1);
t213 = pkin(4) * t363 + t260;
t88 = t294 * t144 - t150 * t290;
t282 = pkin(2) * t358;
t189 = t282 + t302;
t331 = -t189 * t291 + t295 * t246;
t245 = t430 * t358;
t248 = -qJD(2) * pkin(2) + t483;
t328 = m(4) * t248 + (mrSges(4,1) + mrSges(3,3)) * t283 + t479 * qJD(2);
t327 = m(4) * t255 - mrSges(3,3) * t359 + t257 + t371;
t325 = mrSges(5,1) * t291 + mrSges(5,2) * t295;
t324 = Ifges(5,1) * t295 - t395;
t322 = -Ifges(5,2) * t291 + t394;
t62 = pkin(5) * t292 + pkin(10) * t204 + t88;
t63 = pkin(10) * t203 + t89;
t34 = -t289 * t63 + t293 * t62;
t35 = t289 * t62 + t293 * t63;
t158 = mrSges(5,1) * t336 - mrSges(5,3) * t183;
t159 = -mrSges(5,2) * t336 + mrSges(5,3) * t184;
t315 = t295 * t158 + t291 * t159;
t187 = -mrSges(5,2) * t270 + mrSges(5,3) * t235;
t188 = mrSges(5,1) * t270 - mrSges(5,3) * t236;
t313 = t295 * t187 - t291 * t188;
t306 = -qJ(3) * t356 - t355;
t76 = t309 * qJD(2) + (t295 * t334 - t239) * qJD(4) + t331;
t95 = t295 * t189 - t227 * t354 + t291 * t246 + t259 * t353;
t84 = pkin(9) * t303 + t95;
t24 = t144 * t350 - t150 * t351 + t290 * t76 + t294 * t84;
t287 = qJD(2) * qJD(3);
t199 = -qJD(1) * t245 + t287;
t138 = -pkin(4) * t184 + t199;
t25 = -qJD(5) * t89 - t290 * t84 + t294 * t76;
t178 = -pkin(4) * t338 + (-pkin(7) + t340) * t358;
t276 = Ifges(3,4) * t359;
t298 = t130 * mrSges(5,1) + t14 * mrSges(7,1) + t248 * mrSges(4,1) + t50 * mrSges(6,1) + t270 * Ifges(5,3) + t236 * Ifges(5,5) + t235 * Ifges(5,6) + Ifges(3,1) * t283 / 0.2e1 + Ifges(3,5) * t464 + t276 / 0.2e1 + Ifges(4,4) * t465 + (-t292 * Ifges(4,2) - t389) * t466 + t253 * Ifges(7,3) + t105 * Ifges(7,5) + t473 * Ifges(7,6) + t261 * Ifges(6,3) + t164 * Ifges(6,5) + t330 * Ifges(6,6) - t131 * mrSges(5,2) - t15 * mrSges(7,2) - t221 * mrSges(4,3) - t51 * mrSges(6,2);
t268 = Ifges(5,3) * t336;
t247 = pkin(7) * t337 - t287;
t241 = (mrSges(4,2) * t296 - mrSges(4,3) * t292) * qJD(1);
t211 = pkin(4) * t365 + t275 * t289;
t210 = -pkin(4) * t367 + t275 * t293;
t206 = t282 + t306;
t196 = pkin(5) * t311 + t272;
t191 = qJD(1) * t306 + t269;
t170 = -t227 * t291 + t240;
t156 = -pkin(5) * t203 + t213;
t135 = mrSges(6,1) * t261 - mrSges(6,3) * t164;
t134 = -mrSges(6,2) * t261 + mrSges(6,3) * t330;
t124 = pkin(4) * t236 + pkin(5) * t164;
t122 = -mrSges(5,1) * t184 + mrSges(5,2) * t183;
t119 = t183 * Ifges(5,1) + t184 * Ifges(5,4) + Ifges(5,5) * t336;
t118 = t183 * Ifges(5,4) + t184 * Ifges(5,2) + Ifges(5,6) * t336;
t117 = -t310 * t358 + t311 * t439;
t116 = qJD(2) * t304 + t310 * t439;
t107 = -mrSges(6,1) * t330 + mrSges(6,2) * t164;
t96 = -qJD(4) * t171 + t331;
t91 = mrSges(7,1) * t253 - mrSges(7,3) * t105;
t90 = -mrSges(7,2) * t253 + mrSges(7,3) * t473;
t87 = -pkin(5) * t117 + t178;
t71 = -mrSges(6,2) * t336 + mrSges(6,3) * t86;
t70 = mrSges(6,1) * t336 - mrSges(6,3) * t85;
t52 = -pkin(5) * t86 + t138;
t49 = -mrSges(7,1) * t473 + mrSges(7,2) * t105;
t43 = -qJD(6) * t140 - t116 * t289 + t117 * t293;
t42 = qJD(6) * t139 + t116 * t293 + t117 * t289;
t38 = -mrSges(6,1) * t86 + mrSges(6,2) * t85;
t37 = t85 * Ifges(6,1) + t86 * Ifges(6,4) + Ifges(6,5) * t336;
t36 = t85 * Ifges(6,4) + t86 * Ifges(6,2) + Ifges(6,6) * t336;
t27 = -mrSges(7,2) * t336 + mrSges(7,3) * t31;
t26 = mrSges(7,1) * t336 - mrSges(7,3) * t30;
t19 = pkin(10) * t117 + t24;
t18 = pkin(5) * t356 - pkin(10) * t116 + t25;
t17 = t293 * t40 - t377;
t16 = -t289 * t40 - t375;
t8 = -mrSges(7,1) * t31 + mrSges(7,2) * t30;
t7 = t30 * Ifges(7,1) + t31 * Ifges(7,4) + Ifges(7,5) * t336;
t6 = t30 * Ifges(7,4) + t31 * Ifges(7,2) + Ifges(7,6) * t336;
t5 = -qJD(6) * t35 + t18 * t293 - t19 * t289;
t4 = qJD(6) * t34 + t18 * t289 + t19 * t293;
t1 = [((-t251 * mrSges(4,2) + t292 * t342 - 0.2e1 * t429) * qJD(1) + t327 * pkin(7) + t343 * qJD(2) - t476) * t358 + m(7) * (t113 * t87 + t14 * t5 + t15 * t4 + t156 * t52 + t2 * t35 + t3 * t34) + m(6) * (t12 * t89 + t13 * t88 + t138 * t213 + t175 * t178 + t24 * t51 + t25 * t50) + (t266 / 0.2e1 + t267 / 0.2e1 + t268 / 0.2e1 - t191 * mrSges(4,3) + t438 + t440 + t441) * t292 + m(4) * (t191 * t251 + t206 * t221) + (t139 * t2 - t14 * t42 - t140 * t3 + t15 * t43) * mrSges(7,3) + (-t183 * t323 / 0.2e1 - t184 * t321 / 0.2e1 + t199 * t326 + t191 * mrSges(4,2) + t119 * t402 + t118 * t401 - t341 * t247 + (t291 * t75 - t295 * t74) * mrSges(5,3) + (t291 * t420 + t149 * t401 + t270 * (t388 - t390) / 0.2e1 + t324 * t410 + t322 * t411 - t220 * t325 + (t130 * t295 + t131 * t291) * mrSges(5,3)) * qJD(4) + (t298 + t344 * qJD(2) + (-0.2e1 * t428 - t251 * mrSges(4,3) + (-t391 / 0.2e1 - t387 / 0.2e1 - t342) * t296 + Ifges(6,5) * t412 + Ifges(6,6) * t413 + Ifges(7,5) * t421 + Ifges(7,6) * t422 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,2) + 0.3e1 / 0.2e1 * Ifges(3,1) + Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1 + t341 * pkin(7)) * t292) * qJD(1) + t328 * pkin(7)) * qJD(2)) * t296 + t260 * t122 + t206 * t241 - t245 * t174 + t213 * t38 + m(5) * (t130 * t96 + t131 * t95 + t170 * t75 + t171 * t74 + t199 * t260 - t220 * t245) + (-Ifges(6,1) * t204 + Ifges(6,4) * t203) * t432 + (Ifges(7,4) * t140 + Ifges(7,2) * t139) * t435 + (Ifges(7,1) * t140 + Ifges(7,4) * t139) * t436 + (Ifges(6,1) * t116 + Ifges(6,4) * t117) * t416 + (Ifges(6,4) * t116 + Ifges(6,2) * t117) * t418 + t7 * t421 + t6 * t422 + (Ifges(7,1) * t42 + Ifges(7,4) * t43) * t423 + (Ifges(7,4) * t42 + Ifges(7,2) * t43) * t425 + (Ifges(6,5) * t116 + Ifges(6,6) * t117) * t404 + (Ifges(7,5) * t42 + Ifges(7,6) * t43) * t406 + t37 * t412 + t36 * t413 + t117 * t469 + t34 * t26 + t35 * t27 + t43 * t47 / 0.2e1 + t42 * t48 / 0.2e1 + t87 * t49 + t88 * t70 + t89 * t71 + t4 * t90 + t5 * t91 + t113 * (-mrSges(7,1) * t43 + mrSges(7,2) * t42) + t116 * t94 / 0.2e1 + t24 * t134 + t25 * t135 + t52 * (-mrSges(7,1) * t139 + mrSges(7,2) * t140) + t156 * t8 + t170 * t158 + t171 * t159 + t175 * (-mrSges(6,1) * t117 + mrSges(6,2) * t116) + t178 * t107 + t95 * t187 + t96 * t188 + (-t116 * t50 + t117 * t51 + t12 * t203 + t13 * t204) * mrSges(6,3) + t138 * (-mrSges(6,1) * t203 - mrSges(6,2) * t204) + (-Ifges(6,4) * t204 + Ifges(6,2) * t203) * t431; t475 * mrSges(6,3) + (mrSges(6,1) * t362 - mrSges(6,2) * t496) * t175 + t484 * qJD(4) + (mrSges(7,1) * t485 - mrSges(7,2) * t372) * t113 + (Ifges(7,5) * t301 - Ifges(7,6) * t66) * t406 + (-t66 / 0.2e1 - t125 / 0.2e1) * t47 + (Ifges(7,1) * t301 - Ifges(7,4) * t66) * t423 + (Ifges(7,4) * t301 - Ifges(7,2) * t66) * t425 - t443 * mrSges(5,3) + (-t173 / 0.2e1 - t195 / 0.2e1) * t94 + (-Ifges(6,1) * t173 - Ifges(6,4) * t172) * t416 + (-Ifges(6,4) * t173 - Ifges(6,2) * t172) * t418 + (-Ifges(6,5) * t173 - Ifges(6,6) * t172) * t404 + t138 * (mrSges(6,1) * t311 - mrSges(6,2) * t310) + (-Ifges(6,4) * t310 - Ifges(6,2) * t311) * t431 + (-Ifges(6,1) * t310 - Ifges(6,4) * t311) * t432 + (Ifges(7,1) * t445 + Ifges(7,4) * t166) * t436 + (Ifges(7,4) * t445 + Ifges(7,2) * t166) * t435 + t52 * (-mrSges(7,1) * t166 + mrSges(7,2) * t445) + (((-qJ(3) * mrSges(4,1) + t343) * qJD(2) + (t429 + (Ifges(3,4) / 0.2e1 + Ifges(4,6) / 0.2e1) * t292) * qJD(1) + (-t327 + t371) * pkin(7) + t476) * t292 + (-t298 + (Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 - Ifges(4,2) / 0.2e1 + Ifges(3,2) / 0.2e1) * t283 - t276 / 0.2e1 + (t428 - t389 / 0.2e1) * qJD(1) + ((-m(4) * pkin(2) + t479) * qJD(2) - t328) * pkin(7) + (Ifges(6,5) * t408 + Ifges(6,6) * t409 + Ifges(7,5) * t414 + Ifges(7,6) * t415 + t390 / 0.2e1 - t388 / 0.2e1 - pkin(2) * mrSges(4,1) + t344) * qJD(2)) * t296) * qJD(1) + ((-m(5) * t316 + t313) * qJD(4) + m(5) * t443 + t315) * t297 + t444 * t107 + (t199 * qJ(3) - t130 * t151 - t131 * t152 + t220 * t442) * m(5) + m(4) * (-qJ(3) * t247 - qJD(3) * t255) + t493 * mrSges(7,3) + (-t172 / 0.2e1 - t194 / 0.2e1) * t93 + t295 * t119 / 0.2e1 + t272 * t38 - t243 * t174 - t247 * mrSges(4,3) + t6 * t415 + (Ifges(6,1) * t195 + Ifges(6,4) * t194) * t417 + (Ifges(6,4) * t195 + Ifges(6,2) * t194) * t419 + (Ifges(6,5) * t195 + Ifges(6,6) * t194) * t405 + t37 * t408 + t36 * t409 + t7 * t414 + t118 * t402 + t446 * t49 + t449 * t134 + t450 * t135 + (t12 * t177 + t13 * t176 + t138 * t272 + t175 * t444 + t449 * t51 + t450 * t50) * m(6) + t462 * t91 + t463 * t90 + (t113 * t446 + t14 * t462 + t15 * t463 + t196 * t52 + t2 * t81 + t3 * t80) * m(7) + t184 * t322 / 0.2e1 + t183 * t324 / 0.2e1 + t199 * t325 + t80 * t26 + t81 * t27 + qJ(3) * t122 + t360 * qJD(3) + t176 * t70 + t177 * t71 + (-m(4) * t221 - t241) * (-qJ(3) * t359 + t278) + (Ifges(7,1) * t332 + Ifges(7,4) * t125) * t424 - t152 * t187 + (Ifges(7,4) * t332 + Ifges(7,2) * t125) * t426 - t151 * t188 + (Ifges(7,5) * t332 + Ifges(7,6) * t125) * t407 + (t301 / 0.2e1 - t332 / 0.2e1) * t48 + t196 * t8; t445 * t26 - t166 * t27 + t311 * t71 - t310 * t70 - t372 * t91 + t485 * t90 - t496 * t135 + t362 * t134 + t313 * qJD(4) + (t241 + t313) * t283 + (t341 * t359 - t107 - t360 - t49) * qJD(2) - m(4) * (-qJD(2) * t255 - t221 * t283) + t315 + (-qJD(2) * t113 - t493) * m(7) + (-qJD(2) * t175 - t475) * m(6) + (-qJD(2) * t220 - t270 * t316 + t443) * m(5); t164 * t469 + t477 * t90 + (-t113 * t124 + t14 * t478 + t15 * t477 + t2 * t211 + t210 * t3) * m(7) + t478 * t91 - m(6) * (t50 * t57 + t51 * t58) + t498 + t268 + (-Ifges(5,2) * t236 + t149 + t228) * t411 + t438 + (-t236 * t107 + t290 * t71 + t294 * t70 + (t134 * t294 - t135 * t290) * qJD(5) + (t12 * t290 + t13 * t294 - t175 * t236 + t350 * t51 - t351 * t50) * m(6)) * pkin(4) + (t130 * t235 + t131 * t236) * mrSges(5,3) - t220 * (mrSges(5,1) * t236 + mrSges(5,2) * t235) + t211 * t27 + t210 * t26 + t236 * t420 + (Ifges(5,1) * t235 - t396) * t410 + (Ifges(5,5) * t235 - Ifges(5,6) * t236) * t403 - t124 * t49 - t58 * t134 - t57 * t135 - t130 * t187 + t131 * t188; -m(7) * (t14 * t16 + t15 * t17) + (-t164 * t49 + t293 * t26 + t289 * t27 + (-t289 * t91 + t293 * t90) * qJD(6) + (-t113 * t164 - t14 * t349 + t15 * t348 + t2 * t289 + t293 * t3) * m(7)) * pkin(5) + t93 * t416 - t17 * t90 - t16 * t91 - t50 * t134 + t51 * t135 + t498; -t14 * t90 + t15 * t91 + t501;];
tauc  = t1(:);
