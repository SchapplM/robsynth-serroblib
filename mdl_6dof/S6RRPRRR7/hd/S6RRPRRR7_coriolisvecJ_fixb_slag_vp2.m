% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RRPRRR7
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
% Datum: 2019-03-09 14:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RRPRRR7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 13:55:45
% EndTime: 2019-03-09 13:56:12
% DurationCPUTime: 12.48s
% Computational Cost: add. (12868->647), mult. (28897->866), div. (0->0), fcn. (19239->8), ass. (0->307)
t264 = sin(qJ(4));
t268 = cos(qJ(4));
t270 = -pkin(2) - pkin(3);
t212 = -t264 * qJ(3) + t268 * t270;
t165 = t268 * qJD(3) + qJD(4) * t212;
t267 = cos(qJ(5));
t269 = cos(qJ(2));
t336 = qJD(1) * t269;
t265 = sin(qJ(2));
t337 = qJD(1) * t265;
t184 = -t264 * t337 - t268 * t336;
t186 = -t264 * t336 + t268 * t337;
t138 = pkin(4) * t186 - pkin(9) * t184;
t241 = qJ(3) * t336;
t171 = t270 * t337 + t241;
t103 = -t138 + t171;
t248 = pkin(7) * t337;
t205 = pkin(8) * t337 - t248;
t249 = pkin(7) * t336;
t206 = -pkin(8) * t336 + t249;
t141 = t205 * t268 + t206 * t264;
t263 = sin(qJ(5));
t63 = t263 * t103 + t267 * t141;
t469 = t165 * t267 - t63;
t62 = t267 * t103 - t141 * t263;
t468 = -t165 * t263 - t62;
t467 = -qJD(3) + t205;
t193 = t264 * t265 + t268 * t269;
t410 = qJD(2) - qJD(4);
t145 = t410 * t193;
t135 = t145 * qJD(1);
t334 = qJD(2) * t265;
t313 = qJD(1) * t334;
t330 = qJD(4) * t268;
t331 = qJD(4) * t264;
t332 = qJD(2) * t269;
t412 = t264 * t332 + t265 * t330 - t269 * t331;
t136 = qJD(1) * t412 - t268 * t313;
t266 = cos(qJ(6));
t177 = qJD(5) - t184;
t151 = t186 * t267 - t263 * t410;
t188 = -qJD(1) * pkin(1) - pkin(2) * t336 - qJ(3) * t337;
t161 = pkin(3) * t336 - t188;
t102 = -pkin(4) * t184 - pkin(9) * t186 + t161;
t315 = t270 * qJD(2);
t164 = t315 - t467;
t260 = qJD(2) * qJ(3);
t187 = t206 + t260;
t124 = t164 * t264 + t187 * t268;
t111 = -pkin(9) * t410 + t124;
t52 = t267 * t102 - t111 * t263;
t45 = -pkin(10) * t151 + t52;
t35 = pkin(5) * t177 + t45;
t262 = sin(qJ(6));
t150 = -t186 * t263 - t267 * t410;
t53 = t102 * t263 + t111 * t267;
t46 = pkin(10) * t150 + t53;
t356 = t262 * t46;
t14 = t266 * t35 - t356;
t327 = qJD(6) * t266;
t328 = qJD(5) * t267;
t345 = t262 * t263;
t409 = qJD(5) + qJD(6);
t142 = -t266 * t328 - t267 * t327 + t345 * t409;
t355 = t266 * t46;
t15 = t262 * t35 + t355;
t192 = -t266 * t267 + t345;
t195 = t262 * t267 + t263 * t266;
t329 = qJD(5) * t263;
t305 = t265 * t315;
t252 = t265 * qJD(3);
t339 = qJD(1) * t252 + qJD(2) * t241;
t148 = qJD(1) * t305 + t339;
t60 = pkin(4) * t136 - pkin(9) * t135 + t148;
t123 = t164 * t268 - t264 * t187;
t395 = pkin(7) - pkin(8);
t209 = t395 * t334;
t259 = qJD(2) * qJD(3);
t169 = -qJD(1) * t209 + t259;
t228 = t395 * t269;
t210 = qJD(2) * t228;
t289 = qJD(1) * t210;
t68 = qJD(4) * t123 + t268 * t169 + t264 * t289;
t12 = t102 * t328 - t111 * t329 + t263 * t60 + t267 * t68;
t83 = -qJD(5) * t151 - t135 * t263;
t10 = pkin(10) * t83 + t12;
t421 = qJD(6) * t14;
t13 = -qJD(5) * t53 - t263 * t68 + t267 * t60;
t82 = qJD(5) * t150 + t135 * t267;
t9 = pkin(5) * t136 - pkin(10) * t82 + t13;
t2 = t10 * t266 + t262 * t9 + t421;
t222 = Ifges(6,5) * t263 + Ifges(6,6) * t267;
t368 = Ifges(6,4) * t263;
t223 = Ifges(6,2) * t267 + t368;
t367 = Ifges(6,4) * t267;
t224 = Ifges(6,1) * t263 + t367;
t3 = -qJD(6) * t15 - t10 * t262 + t266 * t9;
t303 = mrSges(6,1) * t267 - mrSges(6,2) * t263;
t419 = t192 * t184;
t340 = t419 - t142;
t418 = t195 * t184;
t449 = t409 * t195;
t342 = t418 - t449;
t69 = qJD(4) * t124 + t169 * t264 - t268 * t289;
t36 = -pkin(5) * t83 + t69;
t379 = t195 / 0.2e1;
t380 = -t192 / 0.2e1;
t400 = t83 / 0.2e1;
t401 = t82 / 0.2e1;
t168 = qJD(6) + t177;
t310 = t266 * t150 - t151 * t262;
t92 = t150 * t262 + t151 * t266;
t376 = Ifges(7,4) * t92;
t42 = Ifges(7,2) * t310 + Ifges(7,6) * t168 + t376;
t85 = Ifges(7,4) * t310;
t43 = Ifges(7,1) * t92 + Ifges(7,5) * t168 + t85;
t110 = pkin(4) * t410 - t123;
t294 = Ifges(6,5) * t267 - Ifges(6,6) * t263;
t297 = -Ifges(6,2) * t263 + t367;
t300 = Ifges(6,1) * t267 - t368;
t302 = mrSges(6,1) * t263 + mrSges(6,2) * t267;
t386 = t151 / 0.2e1;
t440 = t300 * t386 + t297 * t150 / 0.2e1 + t294 * t177 / 0.2e1 + t110 * t302;
t31 = qJD(6) * t310 + t262 * t83 + t266 * t82;
t32 = -qJD(6) * t92 - t262 * t82 + t266 * t83;
t6 = t31 * Ifges(7,4) + t32 * Ifges(7,2) + Ifges(7,6) * t136;
t7 = t31 * Ifges(7,1) + t32 * Ifges(7,4) + t136 * Ifges(7,5);
t79 = -pkin(5) * t150 + t110;
t466 = -(t14 * t340 - t15 * t342 + t192 * t2 + t195 * t3) * mrSges(7,3) - (mrSges(7,1) * t342 - mrSges(7,2) * t340) * t79 - (t449 / 0.2e1 - t418 / 0.2e1) * t42 - (t142 / 0.2e1 - t419 / 0.2e1) * t43 + (t222 / 0.2e1 + Ifges(7,5) * t379 + Ifges(7,6) * t380 - Ifges(5,6)) * t136 + t36 * (mrSges(7,1) * t192 + mrSges(7,2) * t195) + (-t303 - mrSges(5,1)) * t69 - t68 * mrSges(5,2) + Ifges(5,5) * t135 + t223 * t400 + t224 * t401 + t379 * t7 + t380 * t6 + t440 * qJD(5);
t213 = t268 * qJ(3) + t264 * t270;
t204 = -pkin(9) + t213;
t371 = pkin(10) - t204;
t311 = qJD(5) * t371;
t348 = t184 * t263;
t325 = pkin(10) * t348;
t465 = t263 * t311 - t325 + t469;
t375 = pkin(10) * t267;
t288 = -pkin(5) * t186 + t184 * t375;
t464 = t267 * t311 - t288 + t468;
t461 = -mrSges(3,1) - mrSges(4,1);
t394 = -pkin(10) - pkin(9);
t316 = qJD(5) * t394;
t71 = t267 * t123 + t263 * t138;
t460 = t263 * t316 + t325 - t71;
t70 = -t123 * t263 + t267 * t138;
t459 = t267 * t316 + t288 - t70;
t447 = -qJD(4) * t213 - t268 * t206 + t264 * t467;
t133 = Ifges(7,3) * t136;
t385 = -t168 / 0.2e1;
t397 = -t92 / 0.2e1;
t399 = -t310 / 0.2e1;
t458 = t133 + (Ifges(7,5) * t310 - Ifges(7,6) * t92) * t385 + (t14 * t310 + t15 * t92) * mrSges(7,3) + (-Ifges(7,2) * t92 + t43 + t85) * t399 - t79 * (mrSges(7,1) * t92 + mrSges(7,2) * t310) + (Ifges(7,1) * t310 - t376) * t397;
t456 = Ifges(5,1) / 0.2e1;
t455 = -Ifges(5,2) / 0.2e1;
t454 = -t336 / 0.2e1;
t162 = t371 * t263;
t163 = t371 * t267;
t105 = t162 * t266 + t163 * t262;
t453 = qJD(6) * t105 + t464 * t262 + t266 * t465;
t106 = t162 * t262 - t163 * t266;
t452 = -qJD(6) * t106 - t262 * t465 + t464 * t266;
t333 = qJD(2) * t268;
t183 = -t263 * t333 + t267 * t337;
t185 = t263 * t337 + t267 * t333;
t275 = t409 * t192;
t451 = -t183 * t266 + t185 * t262 - t195 * t330 + t264 * t275;
t450 = -t183 * t262 - t185 * t266 - t192 * t330 - t264 * t449;
t323 = pkin(5) * t329;
t326 = pkin(5) * t348;
t448 = -t323 + t326 - t447;
t174 = Ifges(5,4) * t184;
t317 = t174 / 0.2e1;
t377 = t267 / 0.2e1;
t378 = -t263 / 0.2e1;
t431 = t186 * t456;
t369 = Ifges(6,4) * t151;
t77 = Ifges(6,2) * t150 + Ifges(6,6) * t177 + t369;
t149 = Ifges(6,4) * t150;
t78 = Ifges(6,1) * t151 + Ifges(6,5) * t177 + t149;
t271 = t161 * mrSges(5,2) - Ifges(5,5) * t410 + t78 * t377 + t77 * t378 + t317 + t431 + t440;
t446 = t271 + (t455 + t456) * t186;
t445 = t3 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t31 + Ifges(7,6) * t32;
t214 = -qJD(2) * pkin(2) + qJD(3) + t248;
t366 = Ifges(4,5) * t269;
t428 = -qJD(2) / 0.2e1;
t429 = -qJD(1) / 0.2e1;
t432 = Ifges(3,4) * t454;
t435 = -Ifges(3,1) / 0.2e1;
t442 = (m(4) * t214 + (mrSges(4,2) + mrSges(3,3)) * t337 + t461 * qJD(2)) * pkin(7) - t188 * mrSges(4,3) - (t265 * Ifges(4,1) - t366) * t429 - t337 * t435 - t432 + t214 * mrSges(4,2) - (Ifges(4,4) + Ifges(3,5)) * t428;
t218 = t249 + t260;
t220 = mrSges(4,2) * t336 + qJD(2) * mrSges(4,3);
t245 = Ifges(4,5) * t337;
t370 = Ifges(3,4) * t265;
t434 = Ifges(4,6) / 0.2e1;
t441 = -(m(4) * t218 - qJD(2) * mrSges(3,2) + mrSges(3,3) * t336 + t220) * pkin(7) + t188 * mrSges(4,1) + qJD(2) * t434 + Ifges(4,3) * t454 + t245 / 0.2e1 + Ifges(3,6) * t428 + (t269 * Ifges(3,2) + t370) * t429 - t218 * mrSges(4,2);
t439 = Ifges(7,1) * t419 + Ifges(7,4) * t418;
t438 = Ifges(7,4) * t419 + Ifges(7,2) * t418;
t437 = Ifges(7,5) * t419 + Ifges(7,6) * t418;
t433 = t42 / 0.2e1;
t430 = t184 * t455;
t225 = t394 * t263;
t227 = t394 * t267;
t152 = t225 * t266 + t227 * t262;
t427 = qJD(6) * t152 + t262 * t459 + t266 * t460;
t154 = t225 * t262 - t227 * t266;
t426 = -qJD(6) * t154 - t262 * t460 + t266 * t459;
t292 = t263 * t53 + t267 * t52;
t286 = t292 * mrSges(6,3);
t215 = -t269 * pkin(2) - t265 * qJ(3) - pkin(1);
t191 = t269 * pkin(3) - t215;
t194 = -t264 * t269 + t265 * t268;
t121 = pkin(4) * t193 - pkin(9) * t194 + t191;
t226 = t395 * t265;
t155 = t226 * t264 + t228 * t268;
t147 = t267 * t155;
t75 = t263 * t121 + t147;
t413 = t268 * t226 - t228 * t264;
t272 = t161 * mrSges(5,1) + t52 * mrSges(6,1) + t14 * mrSges(7,1) - t53 * mrSges(6,2) - t15 * mrSges(7,2) - Ifges(5,4) * t186 + t151 * Ifges(6,5) + t92 * Ifges(7,5) + Ifges(5,6) * t410 + t150 * Ifges(6,6) + t310 * Ifges(7,6) + t177 * Ifges(6,3) + t168 * Ifges(7,3) + t430;
t411 = t317 - t286;
t405 = t13 * mrSges(6,1) - t12 * mrSges(6,2) + Ifges(6,5) * t82 + Ifges(6,6) * t83 + t445;
t403 = t31 / 0.2e1;
t402 = t32 / 0.2e1;
t398 = t310 / 0.2e1;
t396 = t92 / 0.2e1;
t393 = pkin(1) * mrSges(3,1);
t392 = pkin(1) * mrSges(3,2);
t390 = t136 / 0.2e1;
t389 = -t150 / 0.2e1;
t387 = -t151 / 0.2e1;
t384 = t168 / 0.2e1;
t382 = -t177 / 0.2e1;
t374 = t267 * pkin(5);
t363 = t413 * t69;
t354 = -mrSges(5,1) * t410 + mrSges(6,1) * t150 - mrSges(6,2) * t151 - mrSges(5,3) * t186;
t347 = t194 * t263;
t338 = qJ(3) * t332 + t252;
t335 = qJD(2) * t264;
t322 = Ifges(3,5) / 0.2e1 + Ifges(4,4) / 0.2e1;
t321 = -0.3e1 / 0.2e1 * Ifges(4,5) + 0.3e1 / 0.2e1 * Ifges(3,4);
t319 = -Ifges(3,6) / 0.2e1 + t434;
t318 = m(4) * pkin(7) + mrSges(4,2);
t144 = -t265 * t333 + t412;
t159 = t305 + t338;
t67 = pkin(4) * t144 - pkin(9) * t145 + t159;
t96 = qJD(4) * t413 - t209 * t268 + t210 * t264;
t312 = -t263 * t96 + t267 * t67;
t74 = t267 * t121 - t155 * t263;
t33 = t82 * Ifges(6,4) + t83 * Ifges(6,2) + t136 * Ifges(6,6);
t309 = -t33 / 0.2e1 - t12 * mrSges(6,3);
t34 = Ifges(6,1) * t82 + Ifges(6,4) * t83 + Ifges(6,5) * t136;
t308 = t13 * mrSges(6,3) - t34 / 0.2e1;
t307 = t52 * mrSges(6,3) - t78 / 0.2e1;
t306 = t53 * mrSges(6,3) + t77 / 0.2e1;
t203 = pkin(4) - t212;
t299 = Ifges(7,1) * t142 + Ifges(7,4) * t449;
t298 = Ifges(7,1) * t195 - Ifges(7,4) * t192;
t296 = Ifges(7,4) * t142 + Ifges(7,2) * t449;
t295 = Ifges(7,4) * t195 - Ifges(7,2) * t192;
t293 = Ifges(7,5) * t142 + Ifges(7,6) * t449;
t56 = pkin(5) * t193 - t194 * t375 + t74;
t61 = -pkin(10) * t347 + t75;
t29 = -t262 * t61 + t266 * t56;
t30 = t262 * t56 + t266 * t61;
t291 = t263 * t52 - t267 * t53;
t108 = -mrSges(6,2) * t177 + mrSges(6,3) * t150;
t109 = mrSges(6,1) * t177 - mrSges(6,3) * t151;
t157 = mrSges(5,2) * t410 + mrSges(5,3) * t184;
t287 = t108 * t267 - t109 * t263 + t157;
t285 = t145 * t263 + t194 * t328;
t25 = t121 * t328 - t155 * t329 + t263 * t67 + t267 * t96;
t97 = qJD(4) * t155 - t209 * t264 - t268 * t210;
t279 = -qJD(5) * t292 + t12 * t267 - t13 * t263;
t277 = m(6) * t279;
t244 = -pkin(4) - t374;
t211 = -pkin(7) * t313 + t259;
t197 = (-t269 * mrSges(4,1) - mrSges(4,3) * t265) * qJD(1);
t178 = t203 + t374;
t176 = pkin(2) * t334 - t338;
t173 = t192 * t264;
t172 = t195 * t264;
t160 = pkin(2) * t313 - t339;
t137 = -mrSges(5,1) * t184 + mrSges(5,2) * t186;
t134 = Ifges(6,3) * t136;
t132 = t192 * t194;
t131 = t195 * t194;
t122 = pkin(5) * t347 - t413;
t86 = t124 + t326;
t73 = mrSges(7,1) * t168 - mrSges(7,3) * t92;
t72 = -mrSges(7,2) * t168 + mrSges(7,3) * t310;
t55 = -mrSges(6,2) * t136 + mrSges(6,3) * t83;
t54 = mrSges(6,1) * t136 - mrSges(6,3) * t82;
t50 = pkin(5) * t285 + t97;
t48 = -mrSges(7,1) * t310 + mrSges(7,2) * t92;
t44 = -mrSges(6,1) * t83 + mrSges(6,2) * t82;
t38 = -t145 * t195 + t194 * t275;
t37 = -t145 * t192 - t194 * t449;
t26 = -qJD(5) * t75 + t312;
t20 = -mrSges(7,2) * t136 + mrSges(7,3) * t32;
t19 = mrSges(7,1) * t136 - mrSges(7,3) * t31;
t18 = t266 * t45 - t356;
t17 = -t262 * t45 - t355;
t16 = -pkin(10) * t285 + t25;
t11 = -t145 * t375 + pkin(5) * t144 + (-t147 + (pkin(10) * t194 - t121) * t263) * qJD(5) + t312;
t8 = -mrSges(7,1) * t32 + mrSges(7,2) * t31;
t5 = -qJD(6) * t30 + t11 * t266 - t16 * t262;
t4 = qJD(6) * t29 + t11 * t262 + t16 * t266;
t1 = [(-t160 * mrSges(4,3) + (t319 * qJD(2) + (t215 * mrSges(4,1) - t265 * t321 - 0.2e1 * t393) * qJD(1) + t441) * qJD(2)) * t265 + (-t160 * mrSges(4,1) + t318 * t211 + (t322 * qJD(2) + (-0.2e1 * t392 - t215 * mrSges(4,3) + t321 * t269 + (-0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + t318 * pkin(7)) * t265) * qJD(1) + t442) * qJD(2)) * t269 + (Ifges(7,1) * t37 + Ifges(7,4) * t38) * t396 + (Ifges(7,4) * t37 + Ifges(7,2) * t38) * t398 + (Ifges(7,5) * t37 + Ifges(7,6) * t38) * t384 + m(6) * (t110 * t97 + t12 * t75 + t13 * t74 + t25 * t53 + t26 * t52 - t363) + m(5) * (-t123 * t97 + t124 * t96 + t148 * t191 + t155 * t68 + t159 * t161 - t363) + (-Ifges(7,5) * t132 - Ifges(7,6) * t131) * t390 + (-t131 * t2 + t132 * t3 - t14 * t37 + t15 * t38) * mrSges(7,3) + t36 * (mrSges(7,1) * t131 - mrSges(7,2) * t132) + (-Ifges(7,4) * t132 - Ifges(7,2) * t131) * t402 + (-Ifges(7,1) * t132 - Ifges(7,4) * t131) * t403 - t354 * t97 + t176 * t197 + m(7) * (t122 * t36 + t14 * t5 + t15 * t4 + t2 * t30 + t29 * t3 + t50 * t79) + t191 * (mrSges(5,1) * t136 + mrSges(5,2) * t135) + t96 * t157 + t159 * t137 + (-t123 * t145 - t124 * t144 - t135 * t413 - t136 * t155) * mrSges(5,3) - t413 * t44 - t132 * t7 / 0.2e1 - t131 * t6 / 0.2e1 + t122 * t8 + t25 * t108 + t26 * t109 + t79 * (-mrSges(7,1) * t38 + mrSges(7,2) * t37) + t4 * t72 + t5 * t73 + t74 * t54 + t75 * t55 + t50 * t48 + m(4) * (t160 * t215 + t176 * t188) + t37 * t43 / 0.2e1 + t29 * t19 + t30 * t20 + t38 * t433 + (t148 * mrSges(5,2) - Ifges(5,4) * t136 + Ifges(5,1) * t135 + t34 * t377 + t294 * t390 + t297 * t400 + t300 * t401 + t33 * t378 + (mrSges(5,3) + t302) * t69 + (-t12 * t263 - t13 * t267) * mrSges(6,3) + (t222 * t382 + t223 * t389 + t224 * t387 + t110 * t303 + t78 * t378 - t267 * t77 / 0.2e1 + t291 * mrSges(6,3)) * qJD(5)) * t194 + (t148 * mrSges(5,1) - Ifges(5,4) * t135 + t133 / 0.2e1 + t134 / 0.2e1 - t68 * mrSges(5,3) + (Ifges(5,2) + Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t136 + t405) * t193 + (t272 + t430) * t144 + (t431 + t271 + t411) * t145; (-t161 * t171 - t212 * t69 + t213 * t68 + (-t141 + t165) * t124 + t447 * t123) * m(5) + t448 * t48 + ((t432 + (t392 + t366 / 0.2e1) * qJD(1) + (-pkin(2) * mrSges(4,2) + (-m(4) * pkin(2) + t461) * pkin(7) + t322) * qJD(2) - t442) * t269 + (-t245 / 0.2e1 + (t393 + t370 / 0.2e1) * qJD(1) + (-Ifges(4,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + t435) * t336 + (pkin(7) * mrSges(3,2) - qJ(3) * mrSges(4,2) + t319) * qJD(2) - t441) * t265) * qJD(1) + t452 * t73 + (t105 * t3 + t106 * t2 + t14 * t452 + t15 * t453 + t178 * t36 + t448 * t79) * m(7) + t453 * t72 + t437 * t385 + t438 * t399 + t439 * t397 + (t411 + t446) * t184 + t447 * t354 + (-t123 * t184 - t124 * t186 - t135 * t212 - t136 * t213) * mrSges(5,3) + t299 * t396 + t296 * t398 + t293 * t384 + (-t447 * t110 + t203 * t69 + t468 * t52 + t469 * t53) * m(6) + ((-t204 * t109 + t307) * t267 + (-t204 * t108 + t306) * t263) * qJD(5) + (-t204 * t54 + t308) * t263 + (t204 * t55 + t309) * t267 + t203 * t44 + t211 * mrSges(4,3) + qJD(3) * t220 + m(4) * (qJ(3) * t211 + qJD(3) * t218) + (-m(4) * t188 - t197) * (pkin(2) * t337 - t241) + t204 * t277 + t287 * t165 + t178 * t8 - t171 * t137 - t141 * t157 + t105 * t19 + t106 * t20 - t63 * t108 - t62 * t109 - t466 + t272 * t186 - t32 * t295 / 0.2e1 - t31 * t298 / 0.2e1; -qJD(2) * t220 - t185 * t108 - t183 * t109 - t172 * t19 - t173 * t20 + t451 * t73 + t450 * t72 + ((-t137 + t197) * t265 + t318 * t332) * qJD(1) + (-t135 * mrSges(5,3) - qJD(2) * t157 + qJD(4) * t287 - t44 - t8) * t268 + (-t136 * mrSges(5,3) - t263 * t54 + t267 * t55 + (-t108 * t263 - t109 * t267) * qJD(5) + t410 * (-t48 + t354)) * t264 - m(4) * (qJD(2) * t218 - t188 * t337) + (-t172 * t3 - t173 * t2 - t268 * t36 + (t331 - t335) * t79 + t450 * t15 + t451 * t14) * m(7) + ((-qJD(4) * t291 - t69) * t268 + (qJD(4) * t110 + t279) * t264 - t110 * t335 - t183 * t52 - t185 * t53) * m(6) + (-t161 * t337 + t264 * t68 - t268 * t69 - t410 * (-t123 * t264 + t124 * t268)) * m(5); (t293 - t437) * t385 + (t296 - t438) * t399 + (t299 - t439) * t397 + (-t174 / 0.2e1 + t123 * mrSges(5,3) + t286 - t446) * t184 + ((-pkin(9) * t109 - t307) * t267 + (pkin(5) * t48 - pkin(9) * t108 - t306) * t263) * qJD(5) + (-pkin(9) * t54 - t308) * t263 + (pkin(9) * t55 - t309) * t267 + t244 * t8 + t354 * t124 + (t124 * mrSges(5,3) - t272) * t186 + pkin(9) * t277 + t152 * t19 + t154 * t20 - t123 * t157 - t71 * t108 - t70 * t109 - t86 * t48 - pkin(4) * t44 + t295 * t402 + t298 * t403 + (-pkin(4) * t69 - t110 * t124 - t52 * t70 - t53 * t71) * m(6) + t466 + t426 * t73 + t427 * t72 + (t152 * t3 + t154 * t2 + t244 * t36 + (t323 - t86) * t79 + t427 * t15 + t426 * t14) * m(7); (-Ifges(6,2) * t151 + t149 + t78) * t389 - m(7) * (t14 * t17 + t15 * t18) + t77 * t386 + (Ifges(6,1) * t150 - t369) * t387 + (Ifges(6,5) * t150 - Ifges(6,6) * t151) * t382 + t92 * t433 - t110 * (mrSges(6,1) * t151 + mrSges(6,2) * t150) - t52 * t108 + t53 * t109 - t18 * t72 - t17 * t73 + (t150 * t52 + t151 * t53) * mrSges(6,3) + t134 + t405 + (-t151 * t48 + t266 * t19 + t262 * t20 + (-t262 * t73 + t266 * t72) * qJD(6) + (t15 * t327 - t151 * t79 + t266 * t3 + (t2 - t421) * t262) * m(7)) * pkin(5) + t458; -t14 * t72 + t15 * t73 + t42 * t396 + t445 + t458;];
tauc  = t1(:);
