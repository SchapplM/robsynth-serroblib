% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRPRP9
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 17:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRPRP9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:47:06
% EndTime: 2018-11-23 17:47:22
% DurationCPUTime: 16.12s
% Computational Cost: add. (7699->723), mult. (18187->910), div. (0->0), fcn. (11374->6), ass. (0->325)
t427 = Ifges(6,1) + Ifges(7,1);
t437 = -Ifges(6,4) + Ifges(7,5);
t426 = Ifges(7,4) + Ifges(6,5);
t448 = qJD(2) / 0.2e1;
t444 = Ifges(4,1) + Ifges(5,1);
t436 = Ifges(4,5) + Ifges(5,4);
t256 = sin(qJ(5));
t257 = sin(qJ(3));
t259 = cos(qJ(5));
t260 = cos(qJ(3));
t197 = t256 * t257 + t259 * t260;
t416 = qJD(3) - qJD(5);
t129 = t416 * t197;
t261 = cos(qJ(2));
t270 = t197 * t261;
t162 = qJD(1) * t270;
t348 = t129 - t162;
t332 = qJD(5) * t259;
t333 = qJD(5) * t256;
t334 = qJD(3) * t260;
t336 = qJD(3) * t257;
t130 = t256 * t334 + t257 * t332 - t259 * t336 - t260 * t333;
t339 = qJD(1) * t261;
t317 = t260 * t339;
t319 = t257 * t339;
t161 = t256 * t317 - t259 * t319;
t347 = t130 - t161;
t246 = pkin(7) * t339;
t447 = -t257 * qJD(4) - t246 + (t317 - t334) * qJ(4);
t237 = qJD(3) - t339;
t226 = qJD(5) - t237;
t262 = -pkin(3) - pkin(4);
t258 = sin(qJ(2));
t210 = -pkin(2) * t261 - pkin(8) * t258 - pkin(1);
t181 = t210 * qJD(1);
t218 = qJD(2) * pkin(8) + t246;
t126 = t260 * t181 - t257 * t218;
t340 = qJD(1) * t258;
t318 = t260 * t340;
t196 = qJD(2) * t257 + t318;
t92 = pkin(9) * t196 + t126;
t441 = -qJD(4) + t92;
t69 = t237 * t262 - t441;
t222 = t237 * qJ(4);
t127 = t257 * t181 + t260 * t218;
t320 = t257 * t340;
t331 = t260 * qJD(2);
t195 = t320 - t331;
t303 = pkin(9) * t195 + t127;
t82 = t222 + t303;
t19 = -t256 * t82 + t259 * t69;
t422 = qJD(6) - t19;
t17 = -pkin(5) * t226 + t422;
t124 = t195 * t256 + t196 * t259;
t275 = t259 * t195 - t196 * t256;
t217 = -qJD(2) * pkin(2) + pkin(7) * t340;
t106 = t195 * pkin(3) - t196 * qJ(4) + t217;
t86 = -pkin(4) * t195 - t106;
t21 = -pkin(5) * t275 - qJ(6) * t124 + t86;
t116 = Ifges(6,4) * t275;
t361 = Ifges(7,5) * t275;
t424 = t427 * t124 + t426 * t226 + t116 - t361;
t446 = t86 * mrSges(6,2) + t17 * mrSges(7,2) - t19 * mrSges(6,3) - t21 * mrSges(7,3) + t424 / 0.2e1;
t312 = Ifges(3,5) * t448;
t443 = Ifges(4,6) - Ifges(5,6);
t425 = -Ifges(6,6) + Ifges(7,6);
t330 = qJD(1) * qJD(2);
t308 = t258 * t330;
t329 = qJD(2) * qJD(3);
t335 = qJD(3) * t258;
t314 = t257 * t335;
t315 = t261 * t331;
t418 = -t314 + t315;
t149 = qJD(1) * t418 + t260 * t329;
t337 = qJD(2) * t261;
t269 = t257 * t337 + t258 * t334;
t150 = qJD(1) * t269 + t257 * t329;
t42 = qJD(5) * t275 + t149 * t259 + t150 * t256;
t43 = qJD(5) * t124 + t149 * t256 - t259 * t150;
t442 = -t426 * t308 + t427 * t42 + t437 * t43;
t322 = t262 * t257;
t433 = qJD(3) * t322 - t262 * t319 - t447;
t20 = t256 * t69 + t259 * t82;
t18 = qJ(6) * t226 + t20;
t115 = Ifges(7,5) * t124;
t47 = Ifges(7,6) * t226 - Ifges(7,3) * t275 + t115;
t365 = Ifges(6,4) * t124;
t50 = Ifges(6,2) * t275 + Ifges(6,6) * t226 + t365;
t439 = t21 * mrSges(7,1) + t86 * mrSges(6,1) + t47 / 0.2e1 - t50 / 0.2e1 - t18 * mrSges(7,2) - t20 * mrSges(6,3);
t438 = -qJD(2) / 0.2e1;
t378 = Ifges(7,2) + Ifges(6,3);
t352 = t257 * t259;
t273 = t256 * t260 - t352;
t435 = pkin(5) * t347 - qJ(6) * t348 + qJD(6) * t273 + t433;
t434 = t436 * t308 + (-Ifges(4,4) + Ifges(5,5)) * t150 + t444 * t149;
t432 = (-t319 + t336) * pkin(3) + t447;
t362 = Ifges(5,5) * t260;
t367 = Ifges(4,4) * t260;
t431 = t257 * t444 - t362 + t367;
t363 = Ifges(5,5) * t257;
t368 = Ifges(4,4) * t257;
t430 = t260 * t444 + t363 - t368;
t243 = Ifges(3,4) * t339;
t189 = Ifges(5,5) * t196;
t107 = t237 * Ifges(5,6) + t195 * Ifges(5,3) + t189;
t357 = t196 * Ifges(4,4);
t110 = -t195 * Ifges(4,2) + t237 * Ifges(4,6) + t357;
t276 = t126 * t260 + t127 * t257;
t417 = qJD(4) - t126;
t101 = -pkin(3) * t237 + t417;
t103 = t222 + t127;
t278 = t101 * t260 - t103 * t257;
t286 = Ifges(5,3) * t257 + t362;
t290 = -Ifges(4,2) * t257 + t367;
t295 = mrSges(5,1) * t257 - mrSges(5,3) * t260;
t297 = mrSges(4,1) * t257 + mrSges(4,2) * t260;
t359 = Ifges(5,6) * t257;
t360 = Ifges(4,6) * t257;
t364 = Ifges(4,5) * t260;
t366 = Ifges(5,4) * t260;
t383 = t260 / 0.2e1;
t385 = t257 / 0.2e1;
t386 = -t257 / 0.2e1;
t392 = t196 / 0.2e1;
t394 = t195 / 0.2e1;
t395 = -t195 / 0.2e1;
t190 = Ifges(4,4) * t195;
t358 = t195 * Ifges(5,5);
t421 = t196 * t444 + t237 * t436 - t190 + t358;
t428 = t237 / 0.2e1;
t263 = t278 * mrSges(5,2) - t276 * mrSges(4,3) + t106 * t295 + t107 * t385 + t110 * t386 + t217 * t297 + t286 * t394 + t290 * t395 + t430 * t392 + (-t360 + t364 + t359 + t366) * t428 + t421 * t383;
t429 = t263 + Ifges(3,1) * t340 / 0.2e1 + t243 / 0.2e1 + t312;
t63 = pkin(5) * t124 - qJ(6) * t275;
t411 = t42 / 0.2e1;
t409 = t43 / 0.2e1;
t403 = -t275 / 0.2e1;
t404 = t275 / 0.2e1;
t401 = -t124 / 0.2e1;
t390 = -t226 / 0.2e1;
t311 = Ifges(3,6) * t438;
t208 = t259 * qJ(4) + t256 * t262;
t423 = -qJD(5) * t208 + t256 * t441 - t259 * t303;
t351 = t257 * t261;
t238 = pkin(7) * t351;
t254 = t261 * pkin(3);
t381 = pkin(9) * t258;
t113 = pkin(4) * t261 + t238 + t254 + (-t210 - t381) * t260;
t349 = t260 * t261;
t239 = pkin(7) * t349;
t160 = t257 * t210 + t239;
t147 = -qJ(4) * t261 + t160;
t125 = t257 * t381 + t147;
t420 = t256 * t113 + t259 * t125;
t419 = t257 * t436 + t260 * t443;
t301 = pkin(2) * t258 - pkin(8) * t261;
t206 = t301 * qJD(2);
t183 = qJD(1) * t206;
t307 = pkin(7) * t308;
t71 = t181 * t334 + t257 * t183 - t218 * t336 - t260 * t307;
t54 = qJ(4) * t308 + t237 * qJD(4) + t71;
t25 = pkin(9) * t150 + t54;
t382 = pkin(7) * t257;
t321 = -pkin(3) - t382;
t272 = (-pkin(4) + t321) * t258;
t300 = t181 * t336 - t183 * t260 + t218 * t334;
t26 = -pkin(9) * t149 + t272 * t330 + t300;
t4 = -qJD(5) * t20 - t25 * t256 + t259 * t26;
t267 = -pkin(9) * t349 + t272;
t299 = qJD(3) * t239 - t206 * t260 + t210 * t336;
t60 = pkin(9) * t314 + qJD(2) * t267 + t299;
t338 = qJD(2) * t258;
t242 = qJ(4) * t338;
t344 = t257 * t206 + t210 * t334;
t350 = t258 * t260;
t62 = t242 + (-pkin(7) * qJD(2) + pkin(9) * qJD(3)) * t350 + (-qJD(4) + (-pkin(7) * qJD(3) + pkin(9) * qJD(2)) * t257) * t261 + t344;
t9 = -qJD(5) * t420 - t256 * t62 + t259 * t60;
t38 = Ifges(7,6) * t43;
t39 = Ifges(6,6) * t43;
t40 = Ifges(6,5) * t42;
t41 = Ifges(7,4) * t42;
t415 = -t378 * t308 + t38 - t39 + t40 + t41;
t324 = Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t325 = -Ifges(5,6) / 0.2e1 + Ifges(4,6) / 0.2e1;
t326 = Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1;
t327 = -Ifges(6,5) / 0.2e1 - Ifges(7,4) / 0.2e1;
t328 = Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t369 = Ifges(3,4) * t258;
t414 = t324 * t275 + t327 * t124 - t325 * t195 + t328 * t196 - (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t226 + t326 * t237 + t103 * mrSges(5,3) + t126 * mrSges(4,1) + t17 * mrSges(7,1) + t20 * mrSges(6,2) + Ifges(4,6) * t395 + Ifges(5,6) * t394 + t311 - (t261 * Ifges(3,2) + t369) * qJD(1) / 0.2e1 + Ifges(6,6) * t403 + Ifges(7,6) * t404 - t101 * mrSges(5,1) - t127 * mrSges(4,2) - t18 * mrSges(7,3) - t19 * mrSges(6,1) + (Ifges(4,3) + Ifges(5,2)) * t428 + t426 * t401 + t436 * t392 + t378 * t390;
t413 = Ifges(7,5) * t411 - Ifges(7,6) * t308 / 0.2e1 + Ifges(7,3) * t409;
t412 = -t42 * Ifges(6,4) / 0.2e1 + Ifges(6,2) * t409 + Ifges(6,6) * t308 / 0.2e1;
t410 = -t43 / 0.2e1;
t407 = pkin(8) - pkin(9);
t406 = pkin(1) * mrSges(3,1);
t405 = pkin(1) * mrSges(3,2);
t400 = t124 / 0.2e1;
t399 = t149 / 0.2e1;
t398 = -t150 / 0.2e1;
t397 = t150 / 0.2e1;
t393 = -t196 / 0.2e1;
t389 = t226 / 0.2e1;
t388 = -t237 / 0.2e1;
t384 = -t260 / 0.2e1;
t30 = -mrSges(6,1) * t308 - mrSges(6,3) * t42;
t31 = mrSges(7,1) * t308 + t42 * mrSges(7,2);
t377 = t31 - t30;
t29 = -mrSges(7,2) * t43 - mrSges(7,3) * t308;
t32 = mrSges(6,2) * t308 - mrSges(6,3) * t43;
t376 = t32 + t29;
t93 = mrSges(7,2) * t275 + mrSges(7,3) * t226;
t371 = mrSges(6,3) * t275;
t94 = -mrSges(6,2) * t226 + t371;
t375 = t93 + t94;
t370 = mrSges(6,3) * t124;
t95 = mrSges(6,1) * t226 - t370;
t96 = -mrSges(7,1) * t226 + mrSges(7,2) * t124;
t374 = t95 - t96;
t373 = mrSges(4,3) * t195;
t372 = mrSges(4,3) * t196;
t204 = t301 * qJD(1);
t178 = t257 * t204;
t240 = qJ(4) * t340;
t105 = t178 + t240 + (-pkin(7) * t350 + pkin(9) * t351) * qJD(1);
t353 = t204 * t260;
t98 = qJD(1) * t267 - t353;
t45 = t259 * t105 + t256 * t98;
t354 = qJD(2) * mrSges(3,2);
t250 = t257 * qJ(4);
t153 = -mrSges(4,2) * t237 - t373;
t156 = -mrSges(5,2) * t195 + mrSges(5,3) * t237;
t346 = -t153 - t156;
t154 = mrSges(4,1) * t237 - t372;
t155 = -mrSges(5,1) * t237 + mrSges(5,2) * t196;
t345 = -t154 + t155;
t133 = t196 * pkin(3) + t195 * qJ(4);
t342 = qJ(4) * t315 + qJD(4) * t350;
t209 = -t260 * pkin(3) - pkin(2) - t250;
t220 = t407 * t260;
t159 = t210 * t260 - t238;
t188 = t260 * pkin(4) - t209;
t306 = qJD(3) * t220;
t305 = -pkin(7) + t322;
t304 = m(4) * t217 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t195 + mrSges(4,2) * t196 + mrSges(3,3) * t340;
t97 = -pkin(4) * t196 - t133;
t302 = t258 * t321;
t298 = mrSges(4,1) * t260 - mrSges(4,2) * t257;
t296 = mrSges(5,1) * t260 + mrSges(5,3) * t257;
t289 = Ifges(4,2) * t260 + t368;
t285 = -Ifges(5,3) * t260 + t363;
t271 = qJD(2) * t302;
t66 = qJD(1) * t271 + t300;
t280 = t257 * t66 + t260 * t54;
t72 = t257 * t307 - t300;
t279 = -t257 * t72 + t260 * t71;
t44 = -t105 * t256 + t259 * t98;
t207 = -t256 * qJ(4) + t259 * t262;
t58 = t113 * t259 - t125 * t256;
t219 = t407 * t257;
t274 = t259 * t219 - t220 * t256;
t146 = t219 * t256 + t220 * t259;
t152 = -pkin(7) * t318 + t178;
t119 = -mrSges(5,1) * t308 + t149 * mrSges(5,2);
t167 = -qJ(4) * t333 + t259 * qJD(4) + t262 * t332;
t3 = t259 * t25 + t256 * t26 + t69 * t332 - t333 * t82;
t8 = t113 * t332 - t125 * t333 + t256 * t60 + t259 * t62;
t173 = t197 * t258;
t235 = qJ(4) * t350;
t144 = t258 * t305 + t235;
t61 = t150 * pkin(3) - t149 * qJ(4) + qJD(2) * t246 - t196 * qJD(4);
t1 = -qJ(6) * t308 + qJD(6) * t226 + t3;
t2 = pkin(5) * t308 - t4;
t268 = t4 * mrSges(6,1) - t2 * mrSges(7,1) - t3 * mrSges(6,2) + t1 * mrSges(7,3);
t28 = -pkin(4) * t150 - t61;
t89 = (-t258 * t331 - t261 * t336) * pkin(7) + t344;
t266 = -t72 * mrSges(4,1) + t66 * mrSges(5,1) + t71 * mrSges(4,2) - t54 * mrSges(5,3) + t268;
t70 = (t260 * t262 - t250) * t335 + t305 * t337 + t342;
t234 = Ifges(5,2) * t308;
t233 = Ifges(4,3) * t308;
t216 = mrSges(3,3) * t339 - t354;
t205 = t407 * t336;
t203 = pkin(5) - t207;
t202 = -qJ(6) + t208;
t172 = t256 * t350 - t258 * t352;
t166 = -t235 + (pkin(3) * t257 + pkin(7)) * t258;
t163 = -qJD(6) + t167;
t151 = pkin(7) * t320 + t353;
t148 = -t159 + t254;
t142 = Ifges(5,4) * t149;
t141 = Ifges(4,5) * t149;
t140 = Ifges(4,6) * t150;
t139 = Ifges(5,6) * t150;
t134 = mrSges(5,1) * t195 - mrSges(5,3) * t196;
t132 = qJD(1) * t302 - t353;
t131 = t152 + t240;
t120 = -mrSges(4,2) * t308 - mrSges(4,3) * t150;
t118 = mrSges(4,1) * t308 - mrSges(4,3) * t149;
t117 = -mrSges(5,2) * t150 + mrSges(5,3) * t308;
t99 = pkin(5) * t197 + qJ(6) * t273 + t188;
t90 = t338 * t382 - t299;
t88 = pkin(3) * t269 + pkin(7) * t337 + qJ(4) * t314 - t342;
t85 = t271 + t299;
t84 = mrSges(4,1) * t150 + mrSges(4,2) * t149;
t83 = mrSges(5,1) * t150 - mrSges(5,3) * t149;
t81 = -qJD(4) * t261 + t242 + t89;
t80 = qJD(5) * t146 - t205 * t256 - t259 * t306;
t79 = qJD(5) * t274 - t259 * t205 + t256 * t306;
t76 = t149 * Ifges(4,4) - t150 * Ifges(4,2) + Ifges(4,6) * t308;
t75 = t149 * Ifges(5,5) + Ifges(5,6) * t308 + t150 * Ifges(5,3);
t74 = t258 * t273 * t416 + qJD(2) * t270;
t73 = qJD(5) * t173 + t256 * t418 - t259 * t269;
t68 = pkin(5) * t172 - qJ(6) * t173 + t144;
t65 = -mrSges(6,1) * t275 + mrSges(6,2) * t124;
t64 = -mrSges(7,1) * t275 - mrSges(7,3) * t124;
t53 = -pkin(5) * t261 - t58;
t46 = qJ(6) * t261 + t420;
t36 = pkin(5) * t340 - t44;
t35 = -qJ(6) * t340 + t45;
t34 = t256 * t303 + t259 * t92;
t22 = -t63 + t97;
t16 = mrSges(6,1) * t43 + mrSges(6,2) * t42;
t15 = mrSges(7,1) * t43 - mrSges(7,3) * t42;
t10 = pkin(5) * t73 - qJ(6) * t74 - qJD(6) * t173 + t70;
t7 = pkin(5) * t338 - t9;
t6 = -qJ(6) * t338 + qJD(6) * t261 + t8;
t5 = pkin(5) * t43 - qJ(6) * t42 - qJD(6) * t124 + t28;
t11 = [m(6) * (t144 * t28 + t19 * t9 + t20 * t8 + t3 * t420 + t4 * t58 + t70 * t86) + t420 * t32 + (t172 * t437 + t173 * t427) * t411 + (-Ifges(6,2) * t404 + Ifges(7,3) * t403 + t425 * t389 + t437 * t400 + t439) * t73 + t442 * t173 / 0.2e1 + (t324 * t43 - t327 * t42 - t142 / 0.2e1 - t139 / 0.2e1 + t140 / 0.2e1 - t141 / 0.2e1 + t325 * t150 + t41 / 0.2e1 + t40 / 0.2e1 - t39 / 0.2e1 + t38 / 0.2e1 + t266 + (t304 * pkin(7) + (0.3e1 / 0.2e1 * Ifges(3,4) * t261 - 0.2e1 * t405) * qJD(1) + t312 + t429) * qJD(2) - t328 * t149 - t233 / 0.2e1 - t234 / 0.2e1) * t261 + (Ifges(6,4) * t404 + Ifges(7,5) * t403 + t426 * t389 + t427 * t400 + t446) * t74 + (-t1 * t172 + t173 * t2) * mrSges(7,2) + (Ifges(7,5) * t173 + Ifges(7,3) * t172) * t409 + (Ifges(6,4) * t173 - Ifges(6,2) * t172) * t410 + t172 * t412 + t172 * t413 + m(7) * (t1 * t46 + t10 * t21 + t17 * t7 + t18 * t6 + t2 * t53 + t5 * t68) + m(5) * (t101 * t85 + t103 * t81 + t106 * t88 + t147 * t54 + t148 * t66 + t166 * t61) + m(4) * (t126 * t90 + t127 * t89 + t159 * t72 + t160 * t71) + (-t172 * t3 - t173 * t4) * mrSges(6,3) + (t61 * t295 + t75 * t385 + t76 * t386 + pkin(7) * t84 + t290 * t398 + t286 * t397 + (-t257 * t71 - t260 * t72) * mrSges(4,3) + (-t257 * t54 + t260 * t66) * mrSges(5,2) + (t217 * t298 + t289 * t394 + t285 * t395 + t106 * t296 + t110 * t384 + (t126 * t257 - t127 * t260) * mrSges(4,3) + (-t101 * t257 - t103 * t260) * mrSges(5,2) + t431 * t393 + t419 * t388 + t421 * t386) * qJD(3) + ((-0.2e1 * t406 + (-0.3e1 / 0.2e1 * Ifges(3,4) + t364 / 0.2e1 - t360 / 0.2e1 + t366 / 0.2e1 + t359 / 0.2e1) * t258 + t327 * t173 - t324 * t172 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + (m(4) * pkin(7) + t297) * pkin(7) - t326 - t378) * t261) * qJD(1) + t311 - pkin(7) * t216 + t414) * qJD(2) + t430 * t399 + (qJD(3) * t107 + t434) * t383) * t258 + t28 * (mrSges(6,1) * t172 + mrSges(6,2) * t173) + t5 * (mrSges(7,1) * t172 - mrSges(7,3) * t173) + t166 * t83 + t159 * t118 + t160 * t120 + t85 * t155 + t81 * t156 + t89 * t153 + t90 * t154 + t144 * t16 + t147 * t117 + t148 * t119 + t88 * t134 + t6 * t93 + t8 * t94 + t9 * t95 + t7 * t96 + t70 * t65 + t68 * t15 + t10 * t64 + t46 * t29 + t53 * t31 + t58 * t30; (t50 - t47) * (t161 / 0.2e1 - t130 / 0.2e1) - t377 * t274 + (-Ifges(7,5) * t273 + Ifges(7,3) * t197) * t409 + (-Ifges(6,4) * t273 - Ifges(6,2) * t197) * t410 + (-t1 * t197 + t17 * t348 - t18 * t347 - t2 * t273) * mrSges(7,2) + (-t19 * t348 - t197 * t3 - t20 * t347 + t273 * t4) * mrSges(6,3) + t5 * (mrSges(7,1) * t197 + mrSges(7,3) * t273) + t28 * (mrSges(6,1) * t197 - mrSges(6,2) * t273) + t431 * t399 + (-t101 * t132 - t103 * t131 + t106 * t432 + t209 * t61) * m(5) + t432 * t134 + (t274 * t4 + t146 * t3 + t188 * t28 + t433 * t86 + (-t45 + t79) * t20 + (-t44 - t80) * t19) * m(6) + t433 * t65 + (t129 * t426 + t130 * t425) * t389 + (t161 * t425 + t162 * t426) * t390 + (t161 * t437 + t162 * t427) * t401 + (Ifges(6,4) * t129 + Ifges(7,5) * t162 - Ifges(6,2) * t130 + Ifges(7,3) * t161) * t404 + ((t117 + t120) * t260 + (-t118 + t119) * t257 + (-m(4) * t276 + m(5) * t278 + t257 * t346 + t260 * t345) * qJD(3) + m(5) * t280 + m(4) * t279) * pkin(8) - t442 * t273 / 0.2e1 + (t129 * t427 + t130 * t437) * t400 + (t197 * t437 - t273 * t427) * t411 - t61 * t296 - t374 * t80 + t375 * t79 + t376 * t146 + t197 * t412 + t197 * t413 + t285 * t397 + t289 * t398 + t75 * t384 + t76 * t383 + (mrSges(6,1) * t347 + mrSges(6,2) * t348) * t86 + (mrSges(7,1) * t347 - mrSges(7,3) * t348) * t21 - m(4) * (t126 * t151 + t127 * t152) + t263 * qJD(3) + t279 * mrSges(4,3) + t280 * mrSges(5,2) + (Ifges(6,4) * t162 + Ifges(7,5) * t129 - Ifges(6,2) * t161 + Ifges(7,3) * t130) * t403 + t434 * t385 + (t1 * t146 - t274 * t2 + t5 * t99 + t435 * t21 + (-t35 + t79) * t18 + (-t36 + t80) * t17) * m(7) + t435 * t64 + t424 * (-t162 / 0.2e1 + t129 / 0.2e1) + t209 * t83 + t188 * t16 - t132 * t155 - t131 * t156 - t152 * t153 - t151 * t154 + t99 * t15 - t35 * t93 - t45 * t94 - t44 * t95 - t36 * t96 - pkin(2) * t84 + (((t216 + t354) * pkin(7) + (t369 / 0.2e1 + t406) * qJD(1) + t311 + (t197 * t425 - t273 * t426) * t438 + t419 * t448 - t414) * t258 + (((-m(4) * pkin(2) - mrSges(3,1) - t298) * qJD(2) - t304) * pkin(7) + t312 + (t405 + (-Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t258) * qJD(1) - t243 / 0.2e1 - t429) * t261) * qJD(1); -(-Ifges(6,2) * t403 + Ifges(7,3) * t404 + t425 * t390 + t437 * t401 - t439) * t124 + (-t195 * t444 + t107 + t189 - t357) * t393 - t415 + (-t436 * t195 - t196 * t443) * t388 + (-Ifges(6,4) * t403 - Ifges(7,5) * t404 - t426 * t390 - t427 * t401 + t446) * t275 + t142 + t139 - t140 + t141 + (t101 * t195 + t103 * t196) * mrSges(5,2) + (t207 * t4 + t208 * t3 - t86 * t97 + (t167 - t34) * t20 + t423 * t19) * m(6) + t423 * t374 + (t1 * t202 + t2 * t203 - t21 * t22 + (t163 - t34) * t18 - t423 * t17) * m(7) - t375 * t34 + (-Ifges(4,2) * t196 - t190 + t421) * t394 + (Ifges(5,3) * t196 - t358) * t395 + t110 * t392 + (-t345 + t372) * t127 + (t346 - t373) * t126 + (-pkin(3) * t66 + qJ(4) * t54 - t101 * t127 + t103 * t417 - t106 * t133) * m(5) - t266 - t217 * (mrSges(4,1) * t196 - mrSges(4,2) * t195) + t202 * t29 + t203 * t31 + t207 * t30 + t208 * t32 - t106 * (mrSges(5,1) * t196 + mrSges(5,3) * t195) + t167 * t94 + t163 * t93 + qJD(4) * t156 - t133 * t134 + qJ(4) * t117 - pkin(3) * t119 - t97 * t65 - t22 * t64 + t233 + t234; -t237 * t156 + (t134 - t64 - t65) * t196 + (t226 * t375 - t377) * t259 + (-t226 * t374 + t376) * t256 + t119 + (t1 * t256 - t196 * t21 - t2 * t259 + t226 * (t17 * t256 + t18 * t259)) * m(7) + (-t196 * t86 + t256 * t3 + t259 * t4 + t226 * (-t19 * t256 + t20 * t259)) * m(6) + (-t103 * t237 + t106 * t196 + t66) * m(5); t415 + (-pkin(5) * t2 + qJ(6) * t1 - t17 * t20 + t18 * t422 - t21 * t63) * m(7) + (t124 * t18 - t17 * t275) * mrSges(7,2) - t21 * (mrSges(7,1) * t124 - mrSges(7,3) * t275) - t86 * (mrSges(6,1) * t124 + mrSges(6,2) * t275) + (t124 * t425 + t275 * t426) * t390 + (t275 * t427 + t115 - t365 + t47) * t401 + (t370 + t374) * t20 + t50 * t400 + (-Ifges(6,2) * t124 + t116 + t424) * t403 + t268 + qJD(6) * t93 - t63 * t64 - pkin(5) * t31 + qJ(6) * t29 + (Ifges(7,3) * t124 + t361) * t404 + (t371 - t375) * t19; t124 * t64 - t226 * t93 + 0.2e1 * (t2 / 0.2e1 + t21 * t400 + t18 * t390) * m(7) + t31;];
tauc  = t11(:);
