% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
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
% Datum: 2018-11-23 17:05
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:05:15
% EndTime: 2018-11-23 17:05:29
% DurationCPUTime: 14.44s
% Computational Cost: add. (10048->703), mult. (25014->945), div. (0->0), fcn. (17721->8), ass. (0->309)
t277 = sin(pkin(10));
t278 = cos(pkin(10));
t280 = sin(qJ(4));
t283 = cos(qJ(4));
t236 = t277 * t283 + t278 * t280;
t284 = cos(qJ(2));
t295 = t236 * t284;
t200 = qJD(1) * t295;
t219 = t236 * qJD(4);
t341 = t200 - t219;
t348 = t277 * t280;
t235 = -t283 * t278 + t348;
t294 = t235 * t284;
t201 = qJD(1) * t294;
t218 = t235 * qJD(4);
t340 = -t201 + t218;
t369 = pkin(8) + qJ(3);
t246 = t369 * t277;
t247 = t369 * t278;
t332 = qJD(4) * t283;
t335 = qJD(3) * t278;
t336 = qJD(3) * t277;
t138 = -t246 * t332 + t283 * t335 + (-qJD(4) * t247 - t336) * t280;
t281 = sin(qJ(2));
t339 = qJD(1) * t281;
t303 = pkin(2) * t281 - qJ(3) * t284;
t238 = t303 * qJD(1);
t319 = t277 * t339;
t192 = pkin(7) * t319 + t278 * t238;
t344 = t278 * t284;
t299 = pkin(3) * t281 - pkin(8) * t344;
t156 = qJD(1) * t299 + t192;
t220 = t277 * t238;
t345 = t278 * t281;
t346 = t277 * t284;
t296 = -pkin(7) * t345 - pkin(8) * t346;
t176 = qJD(1) * t296 + t220;
t81 = t280 * t156 + t283 * t176;
t71 = qJ(5) * t339 + t81;
t452 = t138 - t71;
t190 = -t280 * t246 + t283 * t247;
t139 = t236 * qJD(3) + qJD(4) * t190;
t80 = t156 * t283 - t280 * t176;
t451 = -t139 - t80;
t244 = -pkin(2) * t284 - qJ(3) * t281 - pkin(1);
t224 = t244 * qJD(1);
t338 = qJD(1) * t284;
t273 = pkin(7) * t338;
t251 = qJD(2) * qJ(3) + t273;
t179 = t278 * t224 - t251 * t277;
t318 = t278 * t339;
t331 = t277 * qJD(2);
t231 = t318 + t331;
t329 = pkin(3) * t338;
t126 = -pkin(8) * t231 + t179 - t329;
t180 = t277 * t224 + t278 * t251;
t297 = -qJD(2) * t278 + t319;
t136 = -pkin(8) * t297 + t180;
t61 = t283 * t126 - t280 * t136;
t418 = qJD(5) - t61;
t434 = Ifges(5,1) + Ifges(6,1);
t433 = -Ifges(5,4) + Ifges(6,5);
t432 = Ifges(6,4) + Ifges(5,5);
t431 = -Ifges(5,6) + Ifges(6,6);
t450 = Ifges(5,3) + Ifges(6,2);
t449 = -pkin(9) * t341 + t452;
t285 = -pkin(4) - pkin(5);
t324 = t285 * t281;
t448 = pkin(9) * t340 - qJD(1) * t324 - t451;
t171 = t280 * t231 + t283 * t297;
t279 = sin(qJ(6));
t282 = cos(qJ(6));
t288 = t283 * t231 - t280 * t297;
t438 = t171 * t279 + t282 * t288;
t85 = t171 * t282 - t279 * t288;
t82 = Ifges(7,4) * t85;
t447 = Ifges(7,2) * t438 - t82;
t446 = -pkin(9) * t288 + t418;
t330 = qJD(1) * qJD(2);
t315 = t281 * t330;
t290 = qJD(2) * t294;
t117 = -qJD(1) * t290 - qJD(4) * t171;
t291 = qJD(2) * t295;
t118 = qJD(1) * t291 + qJD(4) * t288;
t29 = qJD(6) * t85 + t117 * t282 + t118 * t279;
t30 = -qJD(6) * t438 - t117 * t279 + t118 * t282;
t368 = Ifges(7,5) * t29 + Ifges(7,6) * t30;
t266 = qJD(4) - t338;
t259 = qJD(6) - t266;
t384 = -t259 / 0.2e1;
t215 = qJD(2) * t303 - qJD(3) * t281;
t202 = t215 * qJD(1);
t272 = pkin(7) * t339;
t240 = (qJD(3) - t272) * qJD(2);
t154 = t278 * t202 - t240 * t277;
t293 = t299 * qJD(2);
t121 = qJD(1) * t293 + t154;
t155 = t277 * t202 + t278 * t240;
t314 = t284 * t330;
t312 = t277 * t314;
t133 = -pkin(8) * t312 + t155;
t334 = qJD(4) * t280;
t21 = t280 * t121 + t126 * t332 + t283 * t133 - t136 * t334;
t14 = qJ(5) * t315 + t266 * qJD(5) + t21;
t10 = pkin(9) * t118 + t14;
t22 = t121 * t283 - t126 * t334 - t280 * t133 - t136 * t332;
t313 = qJD(2) * t324;
t11 = -pkin(9) * t117 + qJD(1) * t313 - t22;
t40 = t266 * t285 + t446;
t254 = t266 * qJ(5);
t62 = t280 * t126 + t283 * t136;
t46 = pkin(9) * t171 + t62;
t44 = t254 + t46;
t8 = -t279 * t44 + t282 * t40;
t1 = qJD(6) * t8 + t10 * t282 + t11 * t279;
t9 = t279 * t40 + t282 * t44;
t2 = -qJD(6) * t9 - t10 * t279 + t11 * t282;
t415 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t309 = qJD(2) * pkin(2) - qJD(3) - t272;
t191 = pkin(3) * t297 - t309;
t287 = qJ(5) * t288 - t191;
t49 = t171 * t285 + t287;
t445 = (t438 * t9 + t8 * t85) * mrSges(7,3) - Ifges(7,3) * t315 + (Ifges(7,5) * t85 - Ifges(7,6) * t438) * t384 - t49 * (mrSges(7,1) * t438 + mrSges(7,2) * t85) + t368 + t415;
t392 = t118 / 0.2e1;
t444 = Ifges(6,3) * t392;
t225 = t277 * t329 + t273;
t443 = qJ(5) * t340 - qJD(5) * t236 - t225;
t376 = Ifges(7,4) * t438;
t440 = Ifges(7,1) * t85 - t376;
t408 = -Ifges(3,6) / 0.2e1;
t406 = t29 / 0.2e1;
t405 = t30 / 0.2e1;
t39 = Ifges(7,1) * t438 + Ifges(7,5) * t259 + t82;
t437 = t39 / 0.2e1;
t394 = t117 / 0.2e1;
t390 = -t171 / 0.2e1;
t389 = t171 / 0.2e1;
t381 = t266 / 0.2e1;
t387 = t288 / 0.2e1;
t400 = -t438 / 0.2e1;
t436 = -t315 / 0.2e1;
t189 = t283 * t246 + t247 * t280;
t148 = -pkin(9) * t236 + t189;
t149 = pkin(9) * t235 + t190;
t68 = t148 * t279 + t149 * t282;
t430 = -qJD(6) * t68 - t279 * t449 + t282 * t448;
t67 = t148 * t282 - t149 * t279;
t429 = qJD(6) * t67 + t279 * t448 + t282 * t449;
t428 = t117 * t434 + t433 * t118 + t432 * t315;
t169 = Ifges(5,4) * t171;
t356 = Ifges(6,5) * t171;
t427 = t432 * t266 + t288 * t434 - t169 + t356;
t423 = -t285 * t341 - t443;
t241 = -t279 * qJ(5) + t282 * t285;
t422 = qJD(6) * t241 - t279 * t46 + t282 * t446;
t242 = t282 * qJ(5) + t279 * t285;
t421 = -qJD(6) * t242 - t279 * t446 - t282 * t46;
t417 = -pkin(4) * t341 + t443;
t416 = t117 * t432 + t118 * t431 + t315 * t450;
t16 = -pkin(4) * t315 - t22;
t414 = -t22 * mrSges(5,1) + t16 * mrSges(6,1) + t21 * mrSges(5,2) - t14 * mrSges(6,3);
t325 = Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1;
t326 = Ifges(5,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t327 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t57 = -pkin(4) * t266 + t418;
t58 = t254 + t62;
t413 = -t326 * t171 + t325 * t266 + t327 * t288 + t179 * mrSges(4,1) + t58 * mrSges(6,3) + t61 * mrSges(5,1) + t9 * mrSges(7,2) - Ifges(4,6) * t297 / 0.2e1 - Ifges(4,3) * t338 / 0.2e1 + Ifges(4,5) * t231 + qJD(2) * t408 - (Ifges(3,4) * t281 + t284 * Ifges(3,2)) * qJD(1) / 0.2e1 - t259 * Ifges(7,3) - t438 * Ifges(7,5) - t85 * Ifges(7,6) + Ifges(5,6) * t390 + Ifges(6,6) * t389 - pkin(7) * (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t338) - t180 * mrSges(4,2) - t57 * mrSges(6,1) - t62 * mrSges(5,2) - t8 * mrSges(7,1) + t432 * t387 + t450 * t381;
t412 = Ifges(7,4) * t406 + Ifges(7,2) * t405 + Ifges(7,6) * t436;
t411 = Ifges(7,1) * t406 + Ifges(7,4) * t405 + Ifges(7,5) * t436;
t410 = m(4) * pkin(7);
t409 = -Ifges(7,5) / 0.2e1;
t407 = -Ifges(7,6) / 0.2e1;
t404 = Ifges(6,5) * t394 + Ifges(6,6) * t315 / 0.2e1 + t444;
t403 = -t117 * Ifges(5,4) / 0.2e1 + Ifges(5,2) * t392 + Ifges(5,6) * t436;
t402 = -t85 / 0.2e1;
t401 = t85 / 0.2e1;
t399 = t438 / 0.2e1;
t398 = pkin(1) * mrSges(3,1);
t397 = pkin(1) * mrSges(3,2);
t393 = -t118 / 0.2e1;
t388 = -t288 / 0.2e1;
t383 = t259 / 0.2e1;
t382 = -t266 / 0.2e1;
t380 = -t277 / 0.2e1;
t379 = t277 / 0.2e1;
t378 = t278 / 0.2e1;
t367 = mrSges(4,2) * t278;
t366 = mrSges(5,3) * t171;
t365 = mrSges(5,3) * t288;
t364 = Ifges(4,1) * t231;
t363 = Ifges(4,1) * t278;
t362 = Ifges(4,4) * t231;
t361 = Ifges(4,4) * t277;
t360 = Ifges(4,4) * t278;
t359 = Ifges(5,4) * t288;
t357 = Ifges(4,5) * t278;
t355 = Ifges(4,2) * t277;
t354 = Ifges(4,6) * t277;
t134 = t200 * t282 + t201 * t279;
t175 = t235 * t279 + t236 * t282;
t90 = -qJD(6) * t175 + t218 * t279 + t219 * t282;
t352 = t134 - t90;
t135 = t200 * t279 - t201 * t282;
t174 = t235 * t282 - t236 * t279;
t89 = qJD(6) * t174 - t218 * t282 + t219 * t279;
t351 = t135 - t89;
t350 = Ifges(3,5) * qJD(2);
t349 = qJ(5) * t171;
t347 = t277 * t281;
t140 = -mrSges(5,2) * t266 - t366;
t143 = -mrSges(6,2) * t171 + mrSges(6,3) * t266;
t343 = t140 + t143;
t141 = mrSges(5,1) * t266 - t365;
t142 = -mrSges(6,1) * t266 + mrSges(6,2) * t288;
t342 = t141 - t142;
t230 = t278 * t244;
t178 = -pkin(8) * t345 + t230 + (-pkin(7) * t277 - pkin(3)) * t284;
t197 = pkin(7) * t344 + t277 * t244;
t186 = -pkin(8) * t347 + t197;
t102 = t280 * t178 + t283 * t186;
t203 = mrSges(4,1) * t312 + t314 * t367;
t337 = qJD(2) * t281;
t328 = pkin(7) * t337;
t184 = t278 * t215 + t277 * t328;
t226 = (pkin(3) * t331 + pkin(7) * qJD(2)) * t284;
t239 = pkin(3) * t347 + t281 * pkin(7);
t333 = qJD(4) * t281;
t323 = Ifges(4,5) * t338;
t322 = Ifges(4,6) * t338;
t268 = -pkin(3) * t278 - pkin(2);
t55 = t118 * mrSges(5,1) + t117 * mrSges(5,2);
t54 = t118 * mrSges(6,1) - t117 * mrSges(6,3);
t101 = t178 * t283 - t280 * t186;
t311 = -m(4) * t309 - qJD(2) * mrSges(3,1) + mrSges(4,1) * t297 + t231 * mrSges(4,2) + mrSges(3,3) * t339;
t7 = -t30 * mrSges(7,1) + t29 * mrSges(7,2);
t91 = -qJ(5) * t284 + t102;
t211 = t235 * t281;
t307 = -qJ(5) * t211 - t239;
t92 = t284 * pkin(4) - t101;
t306 = mrSges(4,1) * t277 + t367;
t305 = -t361 + t363;
t304 = -t355 + t360;
t64 = pkin(5) * t284 + pkin(9) * t211 + t92;
t210 = t236 * t281;
t66 = pkin(9) * t210 + t91;
t31 = -t279 * t66 + t282 * t64;
t32 = t279 * t64 + t282 * t66;
t69 = -mrSges(7,2) * t259 + mrSges(7,3) * t85;
t70 = mrSges(7,1) * t259 - mrSges(7,3) * t438;
t302 = -t279 * t70 + t282 * t69;
t144 = t210 * t282 + t211 * t279;
t145 = t210 * t279 - t211 * t282;
t298 = qJ(5) * t236 - t268;
t146 = t293 + t184;
t204 = t277 * t215;
t157 = qJD(2) * t296 + t204;
t42 = t146 * t283 - t280 * t157 - t178 * t334 - t186 * t332;
t99 = -mrSges(6,1) * t315 + t117 * mrSges(6,2);
t41 = t280 * t146 + t283 * t157 + t178 * t332 - t186 * t334;
t214 = (pkin(3) * t277 + pkin(7)) * t314;
t150 = -t236 * t333 - t290;
t292 = qJ(5) * t150 - qJD(5) * t211 - t226;
t33 = t118 * pkin(4) - t117 * qJ(5) - qJD(5) * t288 + t214;
t34 = qJ(5) * t337 - qJD(5) * t284 + t41;
t271 = Ifges(3,4) * t338;
t223 = Ifges(3,1) * t339 + t271 + t350;
t213 = (mrSges(4,1) * t281 - mrSges(4,3) * t344) * t330;
t212 = (-mrSges(4,2) * t281 - mrSges(4,3) * t346) * t330;
t199 = -mrSges(4,1) * t338 - mrSges(4,3) * t231;
t198 = mrSges(4,2) * t338 - mrSges(4,3) * t297;
t196 = -pkin(7) * t346 + t230;
t193 = -pkin(7) * t318 + t220;
t188 = (Ifges(4,5) * t281 + t284 * t305) * t330;
t187 = (Ifges(4,6) * t281 + t284 * t304) * t330;
t185 = -t278 * t328 + t204;
t168 = Ifges(6,5) * t288;
t167 = pkin(4) * t235 - t298;
t164 = -Ifges(4,4) * t297 - t323 + t364;
t163 = -Ifges(4,2) * t297 - t322 + t362;
t151 = t332 * t345 - t333 * t348 + t291;
t132 = t235 * t285 + t298;
t127 = pkin(4) * t210 - t307;
t100 = -mrSges(5,2) * t315 - mrSges(5,3) * t118;
t98 = mrSges(5,1) * t315 - mrSges(5,3) * t117;
t97 = -mrSges(6,2) * t118 + mrSges(6,3) * t315;
t96 = mrSges(5,1) * t171 + mrSges(5,2) * t288;
t95 = mrSges(6,1) * t171 - mrSges(6,3) * t288;
t94 = pkin(4) * t288 + t349;
t93 = t210 * t285 + t307;
t77 = -t171 * Ifges(5,2) + t266 * Ifges(5,6) + t359;
t74 = t266 * Ifges(6,6) + t171 * Ifges(6,3) + t168;
t72 = -pkin(4) * t339 - t80;
t65 = t171 * pkin(4) - t287;
t63 = t285 * t288 - t349;
t56 = pkin(4) * t151 - t292;
t48 = -qJD(6) * t145 - t150 * t279 + t151 * t282;
t47 = qJD(6) * t144 + t150 * t282 + t151 * t279;
t43 = -mrSges(7,1) * t85 + mrSges(7,2) * t438;
t38 = Ifges(7,2) * t85 + Ifges(7,6) * t259 + t376;
t36 = t151 * t285 + t292;
t35 = -pkin(4) * t337 - t42;
t26 = pkin(9) * t151 + t34;
t25 = mrSges(7,2) * t315 + mrSges(7,3) * t30;
t24 = -mrSges(7,1) * t315 - mrSges(7,3) * t29;
t23 = -pkin(9) * t150 + t313 - t42;
t15 = pkin(5) * t118 + t33;
t4 = -qJD(6) * t32 + t23 * t282 - t26 * t279;
t3 = qJD(6) * t31 + t23 * t279 + t26 * t282;
t5 = [(Ifges(7,1) * t145 + Ifges(7,4) * t144) * t406 + (mrSges(5,1) * t214 + mrSges(6,1) * t33 - mrSges(6,2) * t14 - mrSges(5,3) * t21 - Ifges(5,2) * t393 + t394 * t433 + t403 + t404 + t444) * t210 + (-t154 * mrSges(4,1) + t155 * mrSges(4,2) + (t163 * t380 + t164 * t378 + t231 * t305 / 0.2e1 + t223 / 0.2e1 + (Ifges(3,5) / 0.2e1 + t304 * t378) * qJD(2) - t309 * t306 + (-t179 * t278 - t180 * t277) * mrSges(4,3) + t311 * pkin(7)) * qJD(2) + t368 / 0.2e1 + Ifges(7,6) * t405 + Ifges(7,5) * t406 - Ifges(6,6) * t392 - Ifges(5,6) * t393 - t432 * t394 + (-0.2e1 * t397 + (-0.3e1 / 0.2e1 * t357 + 0.3e1 / 0.2e1 * t354 + 0.3e1 / 0.2e1 * Ifges(3,4)) * t284 + (Ifges(4,1) * t278 ^ 2 / 0.2e1 - Ifges(7,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) - 0.3e1 / 0.2e1 * Ifges(4,3) + (-0.3e1 / 0.2e1 * t360 + t355) * t277 + (t306 + t410) * pkin(7) - t325) * t281) * t330 + t414 + t415) * t284 + (Ifges(7,4) * t145 + Ifges(7,2) * t144) * t405 + (t427 / 0.2e1 + mrSges(5,2) * t191 + mrSges(6,2) * t57 - mrSges(5,3) * t61 - mrSges(6,3) * t65 + Ifges(5,4) * t390 + Ifges(6,5) * t389 + t381 * t432 + t387 * t434) * t150 + m(4) * (t154 * t196 + t155 * t197 + t179 * t184 + t180 * t185) + t239 * t55 + t226 * t96 + t197 * t212 + t196 * t213 + t185 * t198 + t184 * t199 - t15 * (-mrSges(7,1) * t144 + mrSges(7,2) * t145) + t41 * t140 + t42 * t141 + t35 * t142 + t34 * t143 + t127 * t54 + t101 * t98 + t102 * t100 + t91 * t97 + t92 * t99 + t93 * t7 + t56 * t95 + t4 * t70 + t3 * t69 + t145 * t411 + t144 * t412 + (Ifges(7,1) * t47 + Ifges(7,4) * t48) * t399 + (Ifges(7,4) * t47 + Ifges(7,2) * t48) * t401 + (Ifges(7,5) * t47 + Ifges(7,6) * t48) * t383 + m(7) * (t1 * t32 - t15 * t93 + t2 * t31 + t3 * t9 + t36 * t49 + t4 * t8) + m(6) * (t127 * t33 + t14 * t91 + t16 * t92 + t34 * t58 + t35 * t57 + t56 * t65) + m(5) * (t101 * t22 + t102 * t21 + t191 * t226 + t214 * t239 + t41 * t62 + t42 * t61) + (((Ifges(4,6) * t378 + t408) * qJD(2) + t413) * qJD(2) + pkin(7) * t203 + t187 * t380 + t188 * t378 + (-t154 * t278 - t155 * t277) * mrSges(4,3) + (-0.2e1 * t398 + t145 * t409 + t144 * t407 + (-0.3e1 / 0.2e1 * Ifges(3,4) - t354 + t357 / 0.2e1) * t281 - t327 * t211 - t326 * t210) * t330) * t281 - t416 * t284 / 0.2e1 - t428 * t211 / 0.2e1 + (-t62 * mrSges(5,3) - t58 * mrSges(6,2) + t191 * mrSges(5,1) + t65 * mrSges(6,1) + t74 / 0.2e1 - t77 / 0.2e1 + Ifges(6,3) * t389 - Ifges(5,2) * t390 + t433 * t387 + t431 * t381) * t151 - (mrSges(5,2) * t214 + mrSges(6,2) * t16 - mrSges(5,3) * t22 - mrSges(6,3) * t33 + Ifges(5,4) * t393 + Ifges(6,5) * t392 + t394 * t434) * t211 + t47 * t437 + t31 * t24 + t32 * t25 + t36 * t43 + t48 * t38 / 0.2e1 + t49 * (-mrSges(7,1) * t48 + mrSges(7,2) * t47) + (t1 * t144 - t145 * t2 - t47 * t8 + t48 * t9) * mrSges(7,3); (t187 / 0.2e1 + t155 * mrSges(4,3) + qJD(3) * t198 + qJ(3) * t212) * t278 + (t74 - t77) * (t219 / 0.2e1 - t200 / 0.2e1) + (t188 / 0.2e1 - t154 * mrSges(4,3) - qJD(3) * t199 - qJ(3) * t213) * t277 + ((t350 / 0.2e1 - t223 / 0.2e1 - t271 / 0.2e1 + qJD(1) * t397 + (-t164 / 0.2e1 - t364 / 0.2e1 + t179 * mrSges(4,3) + t309 * mrSges(4,2) + t323 / 0.2e1) * t278 + ((-m(4) * pkin(2) - mrSges(4,1) * t278 - mrSges(3,1)) * qJD(2) - t311) * pkin(7) + (t163 / 0.2e1 + t362 / 0.2e1 + t309 * mrSges(4,1) + t180 * mrSges(4,3) - t322 / 0.2e1 + (pkin(7) * mrSges(4,2) + t363 / 0.2e1 - t361 / 0.2e1) * qJD(2)) * t277) * t284 + ((t398 + (t354 / 0.2e1 + Ifges(3,4) / 0.2e1) * t281) * qJD(1) + (t304 * t379 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t338 + (pkin(7) * mrSges(3,2) + Ifges(4,5) * t379 + t174 * t407 + t175 * t409 - t235 * t326 + t236 * t327 + t408) * qJD(2) - t413) * t281) * qJD(1) + (-Ifges(5,4) * t201 - Ifges(6,5) * t218 - Ifges(5,2) * t200 + Ifges(6,3) * t219) * t389 + (-Ifges(5,4) * t218 - Ifges(6,5) * t201 - Ifges(5,2) * t219 + Ifges(6,3) * t200) * t390 + t427 * (-t218 / 0.2e1 + t201 / 0.2e1) + (t200 * t431 - t201 * t432) * t382 + (t200 * t433 - t201 * t434) * t388 + (-t134 / 0.2e1 + t90 / 0.2e1) * t38 - m(4) * (t179 * t192 + t180 * t193) + (mrSges(7,1) * t352 - mrSges(7,2) * t351) * t49 + (t1 * t174 - t175 * t2 + t351 * t8 - t352 * t9) * mrSges(7,3) + (-t21 * t235 - t22 * t236 + t340 * t61 + t341 * t62) * mrSges(5,3) + (-t14 * t235 + t16 * t236 - t340 * t57 + t341 * t58) * mrSges(6,2) + (-mrSges(6,1) * t341 + mrSges(6,3) * t340) * t65 + (-mrSges(5,1) * t341 - mrSges(5,2) * t340) * t191 - t342 * t139 + t343 * t138 + t268 * t55 + t33 * (mrSges(6,1) * t235 - mrSges(6,3) * t236) + t214 * (mrSges(5,1) * t235 + mrSges(5,2) * t236) - t225 * t96 - pkin(2) * t203 - t193 * t198 - t192 * t199 - t15 * (-mrSges(7,1) * t174 + mrSges(7,2) * t175) + t167 * t54 - t81 * t140 - t80 * t141 - t72 * t142 - t71 * t143 + t132 * t7 + t67 * t24 + t68 * t25 + t175 * t411 + t174 * t412 + (Ifges(7,1) * t89 + Ifges(7,4) * t90) * t399 + (Ifges(7,1) * t135 + Ifges(7,4) * t134) * t400 + (Ifges(7,4) * t89 + Ifges(7,2) * t90) * t401 + (Ifges(7,4) * t135 + Ifges(7,2) * t134) * t402 + t235 * t403 + t235 * t404 + (Ifges(7,4) * t175 + Ifges(7,2) * t174) * t405 + (Ifges(7,1) * t175 + Ifges(7,4) * t174) * t406 + (Ifges(6,5) * t236 + Ifges(6,3) * t235) * t392 + (Ifges(5,4) * t236 - Ifges(5,2) * t235) * t393 + (Ifges(7,5) * t89 + Ifges(7,6) * t90) * t383 + (Ifges(7,5) * t135 + Ifges(7,6) * t134) * t384 + (-t135 / 0.2e1 + t89 / 0.2e1) * t39 + t423 * t43 + t417 * t95 + t428 * t236 / 0.2e1 + t429 * t69 + (t1 * t68 - t132 * t15 + t2 * t67 + t423 * t49 + t429 * t9 + t430 * t8) * m(7) + t430 * t70 + (-t218 * t432 + t219 * t431) * t381 + (-t218 * t434 + t219 * t433) * t387 + (t235 * t433 + t236 * t434) * t394 + (t99 - t98) * t189 + m(4) * (-t179 * t336 + t180 * t335 + (-t154 * t277 + t155 * t278) * qJ(3)) + (t97 + t100) * t190 + (-t189 * t22 + t190 * t21 - t191 * t225 + t214 * t268 + (t138 - t81) * t62 + t451 * t61) * m(5) + (t14 * t190 + t16 * t189 + t167 * t33 + t417 * t65 + t452 * t58 + (t139 - t72) * t57) * m(6); t343 * t171 + t342 * t288 + t314 * t410 - t7 + t231 * t199 + t85 * t69 - t438 * t70 + t297 * t198 - m(4) * (-t179 * t231 - t180 * t297) + t54 + t55 + t203 + (-t438 * t8 + t85 * t9 + t15) * m(7) + (t171 * t58 - t288 * t57 + t33) * m(6) + (t171 * t62 + t288 * t61 + t214) * m(5); (t38 - t440) * t400 - t191 * (mrSges(5,1) * t288 - mrSges(5,2) * t171) - t65 * (mrSges(6,1) * t288 + mrSges(6,3) * t171) + (t171 * t57 + t288 * t58) * mrSges(6,2) + (Ifges(6,3) * t288 - t356) * t390 + t85 * t437 + t447 * t402 - t414 - t445 + t416 + (t342 + t365) * t62 + (-t343 - t366) * t61 + t241 * t24 + t242 * t25 + qJD(5) * t143 + qJ(5) * t97 - pkin(4) * t99 - t94 * t95 - t63 * t43 + t77 * t387 + t421 * t70 + t422 * t69 + (t1 * t242 + t2 * t241 + t421 * t8 + t422 * t9 - t49 * t63) * m(7) + (-pkin(4) * t16 + qJ(5) * t14 + t418 * t58 - t57 * t62 - t65 * t94) * m(6) + (-Ifges(5,2) * t288 - t169 + t427) * t389 + (-t171 * t432 + t288 * t431) * t382 + (-t171 * t434 + t168 - t359 + t74) * t388; t282 * t24 + t279 * t25 + (-t43 + t95) * t288 + t302 * qJD(6) + (-t143 - t302) * t266 + t99 + (t1 * t279 - t288 * t49 + t2 * t282 + t259 * (-t279 * t8 + t282 * t9)) * m(7) + (-t266 * t58 + t288 * t65 + t16) * m(6); t440 * t400 + t38 * t399 - t8 * t69 + t9 * t70 + (t39 - t447) * t402 + t445;];
tauc  = t5(:);
