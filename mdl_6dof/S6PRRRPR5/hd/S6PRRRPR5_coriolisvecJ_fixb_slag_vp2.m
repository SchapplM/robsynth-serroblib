% Calculate vector of centrifugal and Coriolis load on the joints for
% S6PRRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
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
% Datum: 2019-03-08 23:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6PRRRPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 23:23:13
% EndTime: 2019-03-08 23:23:46
% DurationCPUTime: 17.32s
% Computational Cost: add. (13983->735), mult. (38067->1046), div. (0->0), fcn. (30915->14), ass. (0->352)
t266 = sin(pkin(13));
t269 = cos(pkin(13));
t273 = sin(qJ(4));
t277 = cos(qJ(4));
t244 = t266 * t277 + t269 * t273;
t278 = cos(qJ(3));
t267 = sin(pkin(7));
t341 = qJD(2) * t267;
t328 = t278 * t341;
t195 = t244 * t328;
t236 = t244 * qJD(4);
t345 = t195 - t236;
t243 = t266 * t273 - t269 * t277;
t196 = t243 * t328;
t237 = t243 * qJD(4);
t344 = -t196 + t237;
t388 = -qJ(5) - pkin(10);
t322 = qJD(4) * t388;
t234 = qJD(5) * t277 + t273 * t322;
t290 = -qJD(5) * t273 + t277 * t322;
t169 = t269 * t234 + t266 * t290;
t275 = sin(qJ(2));
t268 = sin(pkin(6));
t343 = qJD(1) * t268;
t330 = t275 * t343;
t242 = pkin(9) * t341 + t330;
t274 = sin(qJ(3));
t232 = t274 * t242;
t279 = cos(qJ(2));
t247 = qJD(2) * pkin(2) + t279 * t343;
t270 = cos(pkin(7));
t271 = cos(pkin(6));
t342 = qJD(1) * t271;
t331 = t267 * t342;
t161 = t278 * (t247 * t270 + t331) - t232;
t312 = pkin(3) * t274 - pkin(10) * t278;
t227 = t312 * t341;
t124 = -t161 * t273 + t277 * t227;
t103 = (-qJ(5) * t277 * t278 + pkin(4) * t274) * t341 + t124;
t125 = t277 * t161 + t273 * t227;
t318 = t273 * t328;
t111 = -qJ(5) * t318 + t125;
t52 = t266 * t103 + t269 * t111;
t454 = -t52 + t169;
t253 = qJD(4) - t328;
t396 = t253 / 0.2e1;
t258 = qJD(2) * t270 + qJD(3);
t329 = t274 * t341;
t211 = t258 * t277 - t273 * t329;
t212 = t258 * t273 + t277 * t329;
t321 = t269 * t211 - t212 * t266;
t443 = t321 / 0.2e1;
t297 = t211 * t266 + t269 * t212;
t444 = t297 / 0.2e1;
t414 = Ifges(6,1) * t444 + Ifges(6,4) * t443 + Ifges(6,5) * t396;
t453 = 0.2e1 * t414;
t272 = sin(qJ(6));
t276 = cos(qJ(6));
t352 = t270 * t274;
t428 = t278 * t242 + t247 * t352;
t162 = t274 * t331 + t428;
t146 = pkin(10) * t258 + t162;
t257 = t270 * t342;
t313 = -pkin(3) * t278 - pkin(10) * t274;
t175 = t257 + (qJD(2) * t313 - t247) * t267;
t98 = -t146 * t273 + t277 * t175;
t83 = -qJ(5) * t212 + t98;
t74 = pkin(4) * t253 + t83;
t99 = t146 * t277 + t175 * t273;
t84 = qJ(5) * t211 + t99;
t80 = t269 * t84;
t35 = t266 * t74 + t80;
t29 = pkin(11) * t253 + t35;
t145 = -pkin(3) * t258 - t161;
t119 = -pkin(4) * t211 + qJD(5) + t145;
t49 = -pkin(5) * t321 - pkin(11) * t297 + t119;
t11 = -t272 * t29 + t276 * t49;
t339 = qJD(3) * t278;
t325 = t277 * t339;
t337 = qJD(4) * t277;
t338 = qJD(4) * t273;
t182 = t258 * t337 + (-t274 * t338 + t325) * t341;
t326 = t273 * t339;
t183 = -t258 * t338 + (-t274 * t337 - t326) * t341;
t117 = t182 * t266 - t269 * t183;
t118 = t182 * t269 + t183 * t266;
t348 = t275 * t278;
t349 = t274 * t279;
t295 = t270 * t348 + t349;
t289 = t295 * qJD(2);
t340 = qJD(3) * t267;
t327 = t274 * t340;
t317 = t271 * t327;
t121 = t428 * qJD(3) + (t268 * t289 + t317) * qJD(1);
t96 = -pkin(4) * t183 + t121;
t33 = pkin(5) * t117 - pkin(11) * t118 + t96;
t323 = qJD(2) * t340;
t315 = t274 * t323;
t347 = t278 * t279;
t350 = t274 * t275;
t293 = -t270 * t350 + t347;
t288 = t293 * qJD(2);
t354 = t267 * t278;
t332 = t271 * t354;
t316 = qJD(3) * t332;
t351 = t270 * t278;
t120 = (t247 * t351 - t232) * qJD(3) + (t268 * t288 + t316) * qJD(1);
t292 = t312 * qJD(3);
t189 = (t292 + t330) * t341;
t43 = -qJD(4) * t99 - t120 * t273 + t277 * t189;
t26 = pkin(4) * t315 - qJ(5) * t182 - qJD(5) * t212 + t43;
t42 = t277 * t120 - t146 * t338 + t175 * t337 + t273 * t189;
t30 = qJ(5) * t183 + qJD(5) * t211 + t42;
t8 = t266 * t26 + t269 * t30;
t6 = pkin(11) * t315 + t8;
t1 = qJD(6) * t11 + t272 * t33 + t276 * t6;
t452 = t1 * mrSges(7,2);
t12 = t272 * t49 + t276 * t29;
t2 = -qJD(6) * t12 - t272 * t6 + t276 * t33;
t451 = t2 * mrSges(7,1);
t450 = -pkin(11) * t329 + t454;
t144 = pkin(4) * t318 + t162;
t449 = pkin(4) * t338 - pkin(5) * t345 + t344 * pkin(11) - t144;
t204 = t293 * t343;
t320 = t267 * t330;
t172 = -t204 * t273 + t277 * t320;
t174 = t204 * t277 + t273 * t320;
t241 = pkin(2) * t352 + pkin(9) * t354;
t221 = pkin(10) * t270 + t241;
t222 = (-pkin(2) + t313) * t267;
t164 = t277 * t221 + t273 * t222;
t228 = t267 * t292;
t355 = t267 * t274;
t259 = pkin(9) * t355;
t240 = pkin(2) * t351 - t259;
t229 = t240 * qJD(3);
t105 = -qJD(4) * t164 + t277 * t228 - t229 * t273;
t238 = t270 * t277 - t273 * t355;
t193 = qJD(4) * t238 + t267 * t325;
t239 = t270 * t273 + t277 * t355;
t67 = pkin(4) * t327 - qJ(5) * t193 - qJD(5) * t239 + t105;
t104 = -t221 * t338 + t222 * t337 + t273 * t228 + t277 * t229;
t194 = -qJD(4) * t239 - t267 * t326;
t72 = qJ(5) * t194 + qJD(5) * t238 + t104;
t431 = (-t174 + t72) * t269 + (-t172 + t67) * t266;
t230 = t241 * qJD(3);
t159 = -pkin(4) * t194 + t230;
t203 = t295 * t343;
t448 = t159 - t203;
t151 = qJD(6) - t321;
t373 = Ifges(7,3) * t151;
t130 = t253 * t276 - t272 * t297;
t374 = Ifges(7,6) * t130;
t131 = t253 * t272 + t276 * t297;
t377 = Ifges(7,5) * t131;
t44 = t373 + t374 + t377;
t375 = Ifges(6,6) * t253;
t376 = Ifges(6,2) * t321;
t382 = Ifges(6,4) * t297;
t91 = t375 + t376 + t382;
t447 = -t91 / 0.2e1 + t44 / 0.2e1;
t361 = t266 * t84;
t34 = t269 * t74 - t361;
t446 = -t119 * mrSges(6,2) + t34 * mrSges(6,3);
t445 = -t119 * mrSges(6,1) - t11 * mrSges(7,1) + t12 * mrSges(7,2) + t35 * mrSges(6,3);
t442 = -pkin(11) * t327 - t431;
t133 = t193 * t266 - t269 * t194;
t134 = t193 * t269 + t194 * t266;
t441 = pkin(5) * t133 - pkin(11) * t134 + t448;
t168 = t234 * t266 - t269 * t290;
t51 = t103 * t269 - t111 * t266;
t440 = -t51 - t168;
t439 = t273 / 0.2e1;
t265 = -pkin(4) * t277 - pkin(3);
t186 = pkin(5) * t243 - pkin(11) * t244 + t265;
t254 = t388 * t277;
t324 = t388 * t273;
t200 = -t269 * t254 + t266 * t324;
t129 = t186 * t272 + t200 * t276;
t438 = -qJD(6) * t129 - t272 * t450 + t276 * t449;
t128 = t186 * t276 - t200 * t272;
t437 = qJD(6) * t128 + t272 * t449 + t276 * t450;
t436 = pkin(5) * t329 - t440;
t435 = -t258 * Ifges(4,6) / 0.2e1;
t163 = -t221 * t273 + t277 * t222;
t126 = -pkin(4) * t354 - qJ(5) * t239 + t163;
t137 = qJ(5) * t238 + t164;
t70 = t266 * t126 + t269 * t137;
t66 = -pkin(11) * t354 + t70;
t176 = -t269 * t238 + t239 * t266;
t177 = t238 * t266 + t239 * t269;
t220 = t259 + (-pkin(2) * t278 - pkin(3)) * t270;
t181 = -pkin(4) * t238 + t220;
t89 = pkin(5) * t176 - pkin(11) * t177 + t181;
t31 = -t272 * t66 + t276 * t89;
t434 = qJD(6) * t31 + t272 * t441 - t276 * t442;
t32 = t272 * t89 + t276 * t66;
t433 = -qJD(6) * t32 + t272 * t442 + t276 * t441;
t301 = t11 * t276 + t12 * t272;
t432 = t301 * mrSges(7,3);
t430 = t104 - t174;
t429 = t105 - t172;
t427 = t270 * t347 - t350;
t426 = -t273 * t43 + t277 * t42;
t425 = t1 * t276 - t2 * t272;
t365 = t212 * Ifges(5,4);
t139 = t211 * Ifges(5,2) + t253 * Ifges(5,6) + t365;
t207 = Ifges(5,4) * t211;
t140 = t212 * Ifges(5,1) + t253 * Ifges(5,5) + t207;
t298 = t273 * t99 + t277 * t98;
t384 = Ifges(5,4) * t277;
t385 = Ifges(5,4) * t273;
t393 = t277 / 0.2e1;
t400 = t212 / 0.2e1;
t401 = t211 / 0.2e1;
t424 = t298 * mrSges(5,3) - t145 * (mrSges(5,1) * t273 + mrSges(5,2) * t277) - (-Ifges(5,2) * t273 + t384) * t401 - (Ifges(5,1) * t277 - t385) * t400 - (Ifges(5,5) * t277 - Ifges(5,6) * t273) * t396 + t139 * t439 - t140 * t393;
t206 = -t247 * t267 + t257;
t335 = Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1;
t386 = Ifges(4,4) * t274;
t423 = t335 * t253 + t206 * mrSges(4,1) + t34 * mrSges(6,1) + t98 * mrSges(5,1) + t212 * Ifges(5,5) + t211 * Ifges(5,6) + t435 - (Ifges(4,2) * t278 + t386) * t341 / 0.2e1 + t297 * Ifges(6,5) + t321 * Ifges(6,6) - t162 * mrSges(4,3) - t35 * mrSges(6,2) - t99 * mrSges(5,2) + (Ifges(5,3) + Ifges(6,3)) * t396;
t7 = t26 * t269 - t266 * t30;
t422 = -t43 * mrSges(5,1) - t7 * mrSges(6,1) + t42 * mrSges(5,2) + t8 * mrSges(6,2) - Ifges(5,5) * t182 - Ifges(6,5) * t118 - Ifges(5,6) * t183 + Ifges(6,6) * t117;
t58 = -qJD(6) * t131 - t118 * t272 + t276 * t315;
t55 = Ifges(7,6) * t58;
t57 = qJD(6) * t130 + t118 * t276 + t272 * t315;
t56 = Ifges(7,5) * t57;
t13 = Ifges(7,3) * t117 + t55 + t56;
t421 = t13 / 0.2e1;
t15 = t57 * Ifges(7,1) + t58 * Ifges(7,4) + t117 * Ifges(7,5);
t420 = t15 / 0.2e1;
t381 = Ifges(7,4) * t131;
t45 = Ifges(7,2) * t130 + Ifges(7,6) * t151 + t381;
t418 = -t45 / 0.2e1;
t417 = t57 / 0.2e1;
t416 = t58 / 0.2e1;
t413 = t117 / 0.2e1;
t412 = -t130 / 0.2e1;
t411 = t130 / 0.2e1;
t410 = -t131 / 0.2e1;
t409 = t131 / 0.2e1;
t408 = -t151 / 0.2e1;
t407 = t151 / 0.2e1;
t406 = -t176 / 0.2e1;
t405 = t177 / 0.2e1;
t404 = t182 / 0.2e1;
t403 = t183 / 0.2e1;
t399 = t238 / 0.2e1;
t398 = t239 / 0.2e1;
t397 = -t253 / 0.2e1;
t395 = -t272 / 0.2e1;
t394 = t276 / 0.2e1;
t392 = pkin(4) * t212;
t380 = Ifges(7,4) * t272;
t379 = Ifges(7,4) * t276;
t372 = t117 * Ifges(6,4);
t371 = t118 * Ifges(6,1);
t370 = t118 * Ifges(6,4);
t369 = t120 * mrSges(4,2);
t102 = mrSges(6,1) * t315 - mrSges(6,3) * t118;
t18 = -mrSges(7,1) * t58 + mrSges(7,2) * t57;
t358 = t18 - t102;
t136 = mrSges(6,1) * t253 - mrSges(6,3) * t297;
t73 = -mrSges(7,1) * t130 + mrSges(7,2) * t131;
t357 = t73 - t136;
t191 = -t268 * t427 - t332;
t356 = t121 * t191;
t353 = t268 * t275;
t346 = -mrSges(4,1) * t258 - mrSges(5,1) * t211 + mrSges(5,2) * t212 + mrSges(4,3) * t329;
t97 = -mrSges(6,1) * t321 + mrSges(6,2) * t297;
t334 = t97 + t346;
t59 = t117 * mrSges(6,1) + t118 * mrSges(6,2);
t319 = t341 * t353;
t314 = t451 - t452;
t311 = -t1 * t272 - t2 * t276;
t310 = -mrSges(4,1) * t278 + mrSges(4,2) * t274;
t309 = mrSges(7,1) * t276 - mrSges(7,2) * t272;
t308 = mrSges(7,1) * t272 + mrSges(7,2) * t276;
t307 = Ifges(7,1) * t276 - t380;
t306 = Ifges(7,1) * t272 + t379;
t305 = -Ifges(7,2) * t272 + t379;
t304 = Ifges(7,2) * t276 + t380;
t303 = Ifges(7,5) * t276 - Ifges(7,6) * t272;
t302 = Ifges(7,5) * t272 + Ifges(7,6) * t276;
t300 = t11 * t272 - t12 * t276;
t24 = -t266 * t72 + t269 * t67;
t85 = -mrSges(7,2) * t151 + mrSges(7,3) * t130;
t86 = mrSges(7,1) * t151 - mrSges(7,3) * t131;
t299 = -t272 * t86 + t276 * t85;
t294 = t270 * t349 + t348;
t192 = t268 * t294 + t271 * t355;
t235 = -t267 * t268 * t279 + t270 * t271;
t147 = -t192 * t273 + t235 * t277;
t148 = t192 * t277 + t235 * t273;
t94 = t147 * t266 + t148 * t269;
t62 = t191 * t276 - t272 * t94;
t63 = t191 * t272 + t276 * t94;
t69 = t126 * t269 - t137 * t266;
t149 = -t177 * t272 - t276 * t354;
t296 = -t177 * t276 + t272 * t354;
t256 = Ifges(4,4) * t328;
t286 = Ifges(4,1) * t329 / 0.2e1 + t256 / 0.2e1 + t258 * Ifges(4,5) + t206 * mrSges(4,2) - t161 * mrSges(4,3);
t285 = t374 / 0.2e1 + t377 / 0.2e1 + t373 / 0.2e1 - t376 / 0.2e1 - t382 / 0.2e1 - t375 / 0.2e1 + t447;
t28 = -pkin(5) * t253 - t34;
t127 = Ifges(7,4) * t130;
t46 = Ifges(7,1) * t131 + Ifges(7,5) * t151 + t127;
t283 = t28 * t308 + t303 * t407 + t305 * t411 + t307 * t409 + t394 * t46 + t395 * t45;
t281 = t453 + t283;
t264 = -pkin(4) * t269 - pkin(5);
t252 = Ifges(4,5) * t278 * t323;
t251 = Ifges(5,3) * t315;
t250 = Ifges(6,3) * t315;
t226 = t310 * t341;
t225 = -mrSges(4,2) * t258 + mrSges(4,3) * t328;
t218 = (mrSges(4,1) * t274 + mrSges(4,2) * t278) * t323;
t199 = -t254 * t266 - t269 * t324;
t188 = mrSges(5,1) * t253 - mrSges(5,3) * t212;
t187 = -mrSges(5,2) * t253 + mrSges(5,3) * t211;
t173 = -t196 * t276 + t272 * t329;
t171 = t196 * t272 + t276 * t329;
t158 = -mrSges(5,2) * t315 + mrSges(5,3) * t183;
t157 = mrSges(5,1) * t315 - mrSges(5,3) * t182;
t142 = t316 + (qJD(3) * t427 + t288) * t268;
t141 = t317 + (qJD(3) * t294 + t289) * t268;
t135 = -mrSges(6,2) * t253 + mrSges(6,3) * t321;
t122 = -mrSges(5,1) * t183 + mrSges(5,2) * t182;
t110 = t182 * Ifges(5,1) + t183 * Ifges(5,4) + Ifges(5,5) * t315;
t109 = t182 * Ifges(5,4) + t183 * Ifges(5,2) + Ifges(5,6) * t315;
t106 = -t269 * t172 + t174 * t266;
t101 = -mrSges(6,2) * t315 - mrSges(6,3) * t117;
t93 = -t269 * t147 + t148 * t266;
t82 = pkin(5) * t297 - pkin(11) * t321 + t392;
t79 = qJD(4) * t147 + t142 * t277 + t273 * t319;
t78 = -qJD(4) * t148 - t142 * t273 + t277 * t319;
t77 = qJD(6) * t296 - t134 * t272 + t276 * t327;
t76 = qJD(6) * t149 + t134 * t276 + t272 * t327;
t65 = pkin(5) * t354 - t69;
t54 = Ifges(6,5) * t315 + t371 - t372;
t53 = -t117 * Ifges(6,2) + Ifges(6,6) * t315 + t370;
t41 = -mrSges(7,2) * t117 + mrSges(7,3) * t58;
t40 = mrSges(7,1) * t117 - mrSges(7,3) * t57;
t39 = t269 * t83 - t361;
t38 = t266 * t83 + t80;
t37 = t266 * t78 + t269 * t79;
t36 = t266 * t79 - t269 * t78;
t21 = -pkin(5) * t327 - t24;
t17 = t272 * t82 + t276 * t39;
t16 = -t272 * t39 + t276 * t82;
t14 = t57 * Ifges(7,4) + t58 * Ifges(7,2) + t117 * Ifges(7,6);
t10 = -qJD(6) * t63 + t141 * t276 - t272 * t37;
t9 = qJD(6) * t62 + t141 * t272 + t276 * t37;
t5 = -pkin(5) * t315 - t7;
t3 = [t10 * t86 + t94 * t101 + t37 * t135 + t142 * t225 + t147 * t157 + t148 * t158 + t79 * t187 + t78 * t188 + t235 * t218 + t62 * t40 + t63 * t41 + t9 * t85 + t358 * t93 + t357 * t36 + (-mrSges(3,1) * t275 - mrSges(3,2) * t279) * qJD(2) ^ 2 * t268 + (t59 + t122) * t191 + t334 * t141 + (t226 * t353 + (t191 * t278 - t192 * t274) * qJD(3) * mrSges(4,3)) * t341 + m(5) * (t141 * t145 + t147 * t43 + t148 * t42 + t78 * t98 + t79 * t99 + t356) + m(6) * (t119 * t141 + t191 * t96 - t34 * t36 + t35 * t37 - t7 * t93 + t8 * t94) + m(7) * (t1 * t63 + t10 * t11 + t12 * t9 + t2 * t62 + t28 * t36 + t5 * t93) + m(4) * (t120 * t192 + t356 - t141 * t161 + t142 * t162 + (qJD(1) * t235 + t206) * t319); (Ifges(7,5) * t76 + Ifges(7,6) * t77) * t407 + (-t226 * t330 + t121 * mrSges(4,3) * t274 - pkin(2) * t218 + (t120 * mrSges(4,3) - t250 / 0.2e1 - t251 / 0.2e1 + t422) * t278) * t267 + (Ifges(7,4) * t76 + Ifges(7,2) * t77) * t411 + (t181 * t96 + t69 * t7 + t70 * t8 + t431 * t35 + (t106 + t24) * t34 + t448 * t119) * m(6) + (-Ifges(7,5) * t296 + Ifges(7,6) * t149 + Ifges(7,3) * t176) * t413 + (-Ifges(7,4) * t296 + Ifges(7,2) * t149 + Ifges(7,6) * t176) * t416 + (-Ifges(7,1) * t296 + Ifges(7,4) * t149 + Ifges(7,5) * t176) * t417 + t5 * (-mrSges(7,1) * t149 - mrSges(7,2) * t296) - t296 * t420 - t117 * (Ifges(6,4) * t177 - Ifges(6,2) * t176) / 0.2e1 + (Ifges(7,1) * t76 + Ifges(7,4) * t77) * t409 - m(4) * (-t161 * t203 + t162 * t204 + t206 * t320) + (-t193 * t98 + t194 * t99 + t238 * t42 - t239 * t43) * mrSges(5,3) + t346 * t230 - t334 * t203 + (-Ifges(6,4) * t444 + Ifges(7,5) * t409 - Ifges(6,2) * t443 - Ifges(6,6) * t396 + Ifges(7,6) * t411 + Ifges(7,3) * t407 - t445 + t447) * t133 + ((-m(4) * pkin(2) + t310) * t320 + ((Ifges(4,5) * t270 / 0.2e1 - t240 * mrSges(4,3) + 0.3e1 / 0.2e1 * Ifges(4,4) * t354) * t278 + (-Ifges(4,6) * t270 - t241 * mrSges(4,3) + Ifges(5,5) * t398 + Ifges(5,6) * t399 + Ifges(6,5) * t405 + Ifges(6,6) * t406 - 0.3e1 / 0.2e1 * Ifges(4,4) * t355 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) - t335) * t354) * t274) * qJD(3)) * t341 + (t252 / 0.2e1 - t369 - t121 * mrSges(4,1)) * t270 + m(4) * (t120 * t241 - t121 * t240 - t161 * t230 + t162 * t229) + (-t176 * t8 - t177 * t7) * mrSges(6,3) - t357 * t106 + t118 * (Ifges(6,1) * t177 - Ifges(6,4) * t176) / 0.2e1 + t176 * t451 + t176 * t421 + (Ifges(5,4) * t239 + Ifges(5,2) * t238) * t403 + (Ifges(5,1) * t239 + Ifges(5,4) * t238) * t404 + t54 * t405 + t53 * t406 + (Ifges(5,5) * t193 + Ifges(5,6) * t194) * t396 + t110 * t398 + t109 * t399 + (Ifges(5,1) * t193 + Ifges(5,4) * t194) * t400 + (Ifges(5,4) * t193 + Ifges(5,2) * t194) * t401 + (-t204 + t229) * t225 + (t1 * t149 - t11 * t76 + t12 * t77 + t2 * t296) * mrSges(7,3) + t429 * t188 + t430 * t187 + (t121 * t220 + t163 * t43 + t164 * t42 + t430 * t99 + t429 * t98 + (-t203 + t230) * t145) * m(5) + t431 * t135 - t176 * t452 + (t453 - t446) * t134 + t433 * t86 + t434 * t85 + (t1 * t32 + t2 * t31 + t5 * t65 + (-t106 + t21) * t28 + t434 * t12 + t433 * t11) * m(7) + (t286 * t278 + (t435 + t423) * t274) * t340 + t31 * t40 + t32 * t41 + t65 * t18 + t21 * t73 + t76 * t46 / 0.2e1 + t77 * t45 / 0.2e1 + t28 * (-mrSges(7,1) * t77 + mrSges(7,2) * t76) + t70 * t101 + t69 * t102 + t24 * t136 + t149 * t14 / 0.2e1 + t159 * t97 + t163 * t157 + t164 * t158 + t96 * (mrSges(6,1) * t176 + mrSges(6,2) * t177) + t181 * t59 + t193 * t140 / 0.2e1 + t194 * t139 / 0.2e1 + t145 * (-mrSges(5,1) * t194 + mrSges(5,2) * t193) + t220 * t122 + t121 * (-mrSges(5,1) * t238 + mrSges(5,2) * t239); t440 * t136 - t321 * (-Ifges(6,4) * t196 - Ifges(6,2) * t195) / 0.2e1 - t297 * (-Ifges(6,1) * t196 - Ifges(6,4) * t195) / 0.2e1 + (-Ifges(6,5) * t196 - Ifges(6,6) * t195) * t397 + ((t272 * t237 - t171) * mrSges(7,3) + t345 * mrSges(7,2)) * t12 + ((t276 * t237 + t173) * mrSges(7,3) - t345 * mrSges(7,1)) * t11 + ((m(6) * t119 + t97) * t273 * pkin(4) + (-m(5) * t298 - t273 * t187 - t277 * t188) * pkin(10) - t424) * qJD(4) + ((-t256 / 0.2e1 - t286 + t424) * t278 + ((t258 / 0.2e1 - qJD(3)) * Ifges(4,6) + (t386 / 0.2e1 + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t278) * t341 + (Ifges(5,5) * t273 + Ifges(6,5) * t244 + Ifges(5,6) * t277 - Ifges(6,6) * t243) * qJD(3) / 0.2e1 - t423) * t274) * t341 + t252 - t369 + t285 * t236 + (t34 * t344 + t345 * t35) * mrSges(6,3) + (-mrSges(6,1) * t345 - mrSges(6,2) * t344) * t119 - t346 * t162 + (-t7 * mrSges(6,3) + t5 * t308 + t307 * t417 + t305 * t416 + t303 * t413 + t14 * t395 + t15 * t394 + t54 / 0.2e1 + t96 * mrSges(6,2) - t372 / 0.2e1 + t371 / 0.2e1 + t311 * mrSges(7,3) + (mrSges(7,3) * t300 + t276 * t418 + t28 * t309 + t302 * t408 + t304 * t412 + t306 * t410 + t395 * t46) * qJD(6)) * t244 + (-t8 * mrSges(6,3) - t370 / 0.2e1 + t421 + t56 / 0.2e1 + t55 / 0.2e1 - t53 / 0.2e1 + t96 * mrSges(6,1) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t117 + t314) * t243 + t436 * t73 - t281 * t237 - m(5) * (t124 * t98 + t125 * t99 + t145 * t162) - m(6) * (t119 * t144 + t34 * t51 + t35 * t52) + m(6) * (-t168 * t34 + t169 * t35 - t199 * t7 + t200 * t8 + t265 * t96) + t437 * t85 + (t1 * t129 + t11 * t438 + t12 * t437 + t128 * t2 + t199 * t5 + t28 * t436) * m(7) + t438 * t86 + t358 * t199 + t454 * t135 + (-mrSges(5,1) * t277 + mrSges(5,2) * t273 - mrSges(4,1)) * t121 + t171 * t418 + (Ifges(5,2) * t277 + t385) * t403 + (Ifges(5,1) * t273 + t384) * t404 + (Ifges(7,5) * t173 + Ifges(7,6) * t171 + Ifges(7,3) * t195) * t408 + (Ifges(7,1) * t173 + Ifges(7,4) * t171 + Ifges(7,5) * t195) * t410 + (Ifges(7,4) * t173 + Ifges(7,2) * t171 + Ifges(7,6) * t195) * t412 + t109 * t393 + t110 * t439 + (-t157 * t273 + t158 * t277) * pkin(10) + t265 * t59 + m(5) * (-pkin(3) * t121 + pkin(10) * t426) + t426 * mrSges(5,3) - pkin(3) * t122 + t128 * t40 + t129 * t41 - t144 * t97 + t196 * t414 - t173 * t46 / 0.2e1 - t28 * (-mrSges(7,1) * t171 + mrSges(7,2) * t173) - t125 * t187 - t124 * t188 - t195 * t44 / 0.2e1 + t195 * t91 / 0.2e1 + t200 * t101 - t161 * t225; -(-Ifges(5,2) * t212 + t140 + t207) * t211 / 0.2e1 + (-t281 + t432 + t446) * t321 + (-t285 + t445) * t297 + t251 + t250 + (t101 * t266 + t102 * t269 - t212 * t97) * pkin(4) + (t211 * t98 + t212 * t99) * mrSges(5,3) + (-t11 * t16 - t12 * t17 + t264 * t5 - t28 * t38) * m(7) - t212 * (Ifges(5,1) * t211 - t365) / 0.2e1 - t357 * t38 - t5 * t309 + t302 * t413 + t304 * t416 + t306 * t417 + t272 * t420 + t14 * t394 + (Ifges(5,5) * t211 - Ifges(5,6) * t212) * t397 + t139 * t400 + t264 * t18 + (-t272 * t40 + t276 * t41 + m(7) * t425 + (-m(7) * t301 - t272 * t85 - t276 * t86) * qJD(6)) * (pkin(4) * t266 + pkin(11)) + t425 * mrSges(7,3) - t422 + (t283 - t432) * qJD(6) - t17 * t85 - t16 * t86 + (-t119 * t392 + t34 * t38 - t35 * t39 + (t266 * t8 + t269 * t7) * pkin(4)) * m(6) - t39 * t135 - t98 * t187 + t99 * t188 - t145 * (mrSges(5,1) * t212 + mrSges(5,2) * t211); t272 * t41 + t276 * t40 - t357 * t297 + t299 * qJD(6) + (-t135 - t299) * t321 + t59 + (-t151 * t300 - t297 * t28 - t311) * m(7) + (t297 * t34 - t321 * t35 + t96) * m(6); -t28 * (mrSges(7,1) * t131 + mrSges(7,2) * t130) + (Ifges(7,1) * t130 - t381) * t410 + t45 * t409 + (Ifges(7,5) * t130 - Ifges(7,6) * t131) * t408 - t11 * t85 + t12 * t86 + (t11 * t130 + t12 * t131) * mrSges(7,3) + t314 + t13 + (-Ifges(7,2) * t131 + t127 + t46) * t412;];
tauc  = t3(:);
