% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRPRPR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
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
% Datum: 2018-11-23 17:08
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRPRPR12_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPR12_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 17:08:18
% EndTime: 2018-11-23 17:08:35
% DurationCPUTime: 16.84s
% Computational Cost: add. (13429->762), mult. (34625->1047), div. (0->0), fcn. (25815->10), ass. (0->359)
t275 = sin(qJ(2));
t271 = sin(pkin(6));
t351 = qJD(1) * t271;
t329 = t275 * t351;
t248 = pkin(2) * t329;
t278 = cos(qJ(2));
t299 = pkin(9) * t275 - qJ(3) * t278;
t196 = t299 * t351 + t248;
t328 = t278 * t351;
t272 = cos(pkin(6));
t350 = qJD(1) * t272;
t338 = pkin(1) * t350;
t223 = pkin(8) * t328 + t275 * t338;
t198 = pkin(3) * t328 + t223;
t274 = sin(qJ(4));
t277 = cos(qJ(4));
t120 = -t196 * t274 + t198 * t277;
t345 = qJD(4) * t274;
t279 = -pkin(2) - pkin(9);
t356 = qJ(5) - t279;
t474 = -(-qJ(5) * t274 * t275 + pkin(4) * t278) * t351 - t120 - t277 * qJD(5) + t345 * t356;
t121 = t196 * t277 + t198 * t274;
t318 = t277 * t329;
t321 = t277 * t356;
t473 = qJ(5) * t318 + qJD(4) * t321 + qJD(5) * t274 + t121;
t270 = sin(pkin(11));
t359 = t270 * t274;
t361 = cos(pkin(11));
t187 = t318 * t361 - t329 * t359;
t322 = t361 * t277;
t226 = -qJD(4) * t322 + t270 * t345;
t355 = t187 - t226;
t234 = t270 * t277 + t274 * t361;
t188 = t234 * t329;
t227 = t234 * qJD(4);
t352 = t227 + t188;
t253 = t278 * t338;
t472 = qJD(3) - t253;
t462 = t270 * t474 - t361 * t473;
t344 = qJD(4) * t277;
t427 = pkin(3) + pkin(8);
t461 = pkin(4) * t344 - (-pkin(4) * t277 - t427) * t329 + t472;
t471 = -pkin(10) * t328 + t462;
t470 = pkin(5) * t355 + pkin(10) * t352 + t461;
t464 = t270 * t473 + t361 * t474;
t273 = sin(qJ(6));
t276 = cos(qJ(6));
t244 = qJD(4) + t329;
t257 = qJD(2) + t350;
t208 = t257 * t277 - t274 * t328;
t138 = t257 * t279 + t329 * t427 + t472;
t360 = qJ(3) * t275;
t190 = (t278 * t279 - pkin(1) - t360) * t271;
t169 = qJD(1) * t190;
t93 = t138 * t277 - t169 * t274;
t82 = -qJ(5) * t208 + t93;
t76 = pkin(4) * t244 + t82;
t207 = -t257 * t274 - t277 * t328;
t94 = t138 * t274 + t169 * t277;
t83 = qJ(5) * t207 + t94;
t80 = t361 * t83;
t33 = t270 * t76 + t80;
t28 = pkin(10) * t244 + t33;
t246 = t257 * qJ(3);
t161 = t246 + t198;
t122 = -pkin(4) * t207 + qJD(5) + t161;
t289 = t207 * t270 + t208 * t361;
t320 = t207 * t361 - t208 * t270;
t62 = -pkin(5) * t320 - pkin(10) * t289 + t122;
t9 = -t273 * t28 + t276 * t62;
t469 = t9 * mrSges(7,1);
t468 = mrSges(4,1) + mrSges(3,3);
t467 = mrSges(4,2) - mrSges(3,1);
t233 = -t322 + t359;
t265 = pkin(4) * t274 + qJ(3);
t160 = pkin(5) * t234 + pkin(10) * t233 + t265;
t236 = t356 * t274;
t177 = -t236 * t361 - t270 * t321;
t107 = t160 * t276 - t177 * t273;
t466 = qJD(6) * t107 + t273 * t470 + t276 * t471;
t108 = t160 * t273 + t177 * t276;
t465 = -qJD(6) * t108 - t273 * t471 + t276 * t470;
t391 = mrSges(7,3) * t273;
t463 = pkin(5) * t328 - t464;
t366 = t270 * t83;
t32 = t361 * t76 - t366;
t27 = -pkin(5) * t244 - t32;
t310 = mrSges(7,1) * t273 + mrSges(7,2) * t276;
t390 = mrSges(7,3) * t276;
t460 = t27 * t310 - t390 * t9;
t459 = t32 * t352 - t33 * t355;
t222 = pkin(8) * t329 - t253;
t458 = -qJD(3) - t222;
t302 = Ifges(7,5) * t276 - Ifges(7,6) * t273;
t385 = Ifges(7,4) * t276;
t305 = -Ifges(7,2) * t273 + t385;
t386 = Ifges(7,4) * t273;
t308 = Ifges(7,1) * t276 - t386;
t401 = -t276 / 0.2e1;
t403 = t273 / 0.2e1;
t111 = t244 * t276 - t273 * t289;
t128 = qJD(6) - t320;
t112 = t244 * t273 + t276 * t289;
t387 = Ifges(7,4) * t112;
t41 = Ifges(7,2) * t111 + Ifges(7,6) * t128 + t387;
t418 = -t128 / 0.2e1;
t110 = Ifges(7,4) * t111;
t42 = Ifges(7,1) * t112 + Ifges(7,5) * t128 + t110;
t420 = -t112 / 0.2e1;
t422 = -t111 / 0.2e1;
t368 = t244 * Ifges(6,5);
t374 = t289 * Ifges(6,1);
t377 = t320 * Ifges(6,4);
t79 = t368 + t374 + t377;
t428 = -t79 / 0.2e1;
t457 = t302 * t418 + t305 * t422 + t308 * t420 + t401 * t42 + t403 * t41 + t428 - t460;
t10 = t273 * t62 + t276 * t28;
t456 = -t10 * mrSges(7,3) - t41 / 0.2e1;
t405 = -t244 / 0.2e1;
t415 = -t289 / 0.2e1;
t416 = -t320 / 0.2e1;
t455 = Ifges(6,1) * t415 + Ifges(6,4) * t416 + Ifges(6,5) * t405;
t454 = Ifges(6,4) * t415 - Ifges(7,5) * t420 + Ifges(6,2) * t416 + Ifges(6,6) * t405 - Ifges(7,6) * t422 - Ifges(7,3) * t418;
t453 = t244 / 0.2e1;
t452 = -t257 / 0.2e1;
t451 = t257 / 0.2e1;
t450 = -t351 / 0.2e1;
t118 = mrSges(6,1) * t244 - mrSges(6,3) * t289;
t66 = -mrSges(7,1) * t111 + mrSges(7,2) * t112;
t447 = t66 - t118;
t358 = t271 * t275;
t258 = pkin(8) * t358;
t398 = pkin(1) * t278;
t331 = -pkin(2) - t398;
t172 = pkin(3) * t358 + t258 + (-pkin(9) + t331) * t272;
t116 = t172 * t274 + t190 * t277;
t349 = qJD(2) * t271;
t323 = qJD(1) * t349;
t317 = t275 * t323;
t242 = pkin(2) * t317;
t346 = qJD(3) * t275;
t283 = (qJD(2) * t299 - t346) * t271;
t148 = qJD(1) * t283 + t242;
t261 = t272 * t275 * pkin(1);
t357 = t271 * t278;
t199 = (t357 * t427 + t261) * qJD(2);
t173 = qJD(1) * t199;
t60 = t138 * t344 + t148 * t277 - t169 * t345 + t173 * t274;
t61 = -qJD(4) * t94 - t148 * t274 + t173 * t277;
t446 = t274 * t60 + t277 * t61;
t348 = qJD(2) * t275;
t325 = t274 * t348;
t343 = qJD(4) * t278;
t157 = -t257 * t345 + (-t277 * t343 + t325) * t351;
t324 = t277 * t348;
t158 = -t257 * t344 + (t274 * t343 + t324) * t351;
t102 = t157 * t270 - t158 * t361;
t103 = t157 * t361 + t158 * t270;
t316 = t278 * t323;
t58 = qJD(6) * t111 + t103 * t276 + t273 * t316;
t25 = mrSges(7,1) * t102 - mrSges(7,3) * t58;
t59 = -qJD(6) * t112 - t103 * t273 + t276 * t316;
t26 = -mrSges(7,2) * t102 + mrSges(7,3) * t59;
t445 = -t25 * t273 + t26 * t276;
t195 = -t246 - t223;
t300 = -pkin(2) * t278 - t360;
t212 = (-pkin(1) + t300) * t271;
t203 = qJD(1) * t212;
t444 = -(Ifges(3,6) / 0.2e1 - Ifges(4,5) / 0.2e1) * t257 + t195 * mrSges(4,1) + Ifges(3,6) * t452 + (Ifges(3,4) * t275 + Ifges(3,2) * t278) * t450 + Ifges(4,5) * t451 + (-Ifges(4,6) * t275 - Ifges(4,3) * t278) * t351 / 0.2e1 - t203 * mrSges(4,2) - t223 * mrSges(3,3);
t378 = t128 * Ifges(7,3);
t379 = t112 * Ifges(7,5);
t380 = t111 * Ifges(7,6);
t40 = t378 + t379 + t380;
t367 = t244 * Ifges(6,6);
t373 = t289 * Ifges(6,4);
t376 = t320 * Ifges(6,2);
t78 = t367 + t373 + t376;
t443 = -t469 + t78 / 0.2e1 - t40 / 0.2e1;
t31 = pkin(4) * t316 - qJ(5) * t157 - qJD(5) * t208 + t61;
t35 = qJ(5) * t158 + qJD(5) * t207 + t60;
t7 = -t270 * t35 + t31 * t361;
t8 = t270 * t31 + t35 * t361;
t442 = t61 * mrSges(5,1) + t7 * mrSges(6,1) - t60 * mrSges(5,2) - t8 * mrSges(6,2) + Ifges(5,5) * t157 + Ifges(6,5) * t103 + Ifges(5,6) * t158 - Ifges(6,6) * t102;
t186 = -pkin(2) * t257 - t458;
t247 = Ifges(3,4) * t328;
t334 = Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t384 = Ifges(4,6) * t278;
t438 = Ifges(3,5) / 0.2e1;
t441 = t334 * t244 + (t438 - Ifges(4,4) / 0.2e1) * t257 + t186 * mrSges(4,1) + t222 * mrSges(3,3) + t32 * mrSges(6,1) + t93 * mrSges(5,1) + t208 * Ifges(5,5) + t207 * Ifges(5,6) + Ifges(3,1) * t329 / 0.2e1 + Ifges(3,5) * t451 + t247 / 0.2e1 + Ifges(4,4) * t452 + (-Ifges(4,2) * t275 - t384) * t450 + t289 * Ifges(6,5) + t320 * Ifges(6,6) - t203 * mrSges(4,3) - t33 * mrSges(6,2) - t94 * mrSges(5,2) + (Ifges(5,3) + Ifges(6,3)) * t453;
t370 = t208 * Ifges(5,4);
t124 = t207 * Ifges(5,2) + Ifges(5,6) * t244 + t370;
t204 = Ifges(5,4) * t207;
t125 = Ifges(5,1) * t208 + Ifges(5,5) * t244 + t204;
t296 = t274 * t93 - t277 * t94;
t388 = Ifges(5,4) * t277;
t389 = Ifges(5,4) * t274;
t402 = -t274 / 0.2e1;
t409 = -t208 / 0.2e1;
t410 = -t207 / 0.2e1;
t439 = t161 * (mrSges(5,1) * t277 - mrSges(5,2) * t274) + (Ifges(5,5) * t274 + Ifges(5,6) * t277) * t405 + (Ifges(5,2) * t277 + t389) * t410 + (Ifges(5,1) * t274 + t388) * t409 + t296 * mrSges(5,3) + t125 * t402 - t277 * t124 / 0.2e1;
t56 = Ifges(7,6) * t59;
t57 = Ifges(7,5) * t58;
t11 = Ifges(7,3) * t102 + t56 + t57;
t437 = t11 / 0.2e1;
t435 = t40 / 0.2e1;
t433 = t42 / 0.2e1;
t432 = t58 / 0.2e1;
t431 = t59 / 0.2e1;
t430 = -t78 / 0.2e1;
t426 = pkin(1) * mrSges(3,1);
t425 = pkin(1) * mrSges(3,2);
t424 = -t102 / 0.2e1;
t423 = t102 / 0.2e1;
t421 = t111 / 0.2e1;
t419 = t112 / 0.2e1;
t417 = t128 / 0.2e1;
t228 = -t272 * t274 - t277 * t357;
t291 = -t272 * t277 + t274 * t357;
t150 = -t228 * t361 - t270 * t291;
t414 = -t150 / 0.2e1;
t151 = t228 * t270 - t291 * t361;
t413 = t151 / 0.2e1;
t412 = t157 / 0.2e1;
t411 = t158 / 0.2e1;
t408 = t208 / 0.2e1;
t407 = t228 / 0.2e1;
t406 = -t291 / 0.2e1;
t400 = t276 / 0.2e1;
t397 = pkin(4) * t208;
t396 = pkin(4) * t270;
t20 = -mrSges(7,1) * t59 + mrSges(7,2) * t58;
t88 = mrSges(6,1) * t316 - mrSges(6,3) * t103;
t394 = t20 - t88;
t184 = qJD(4) * t228 + t271 * t325;
t347 = qJD(2) * t278;
t326 = t271 * t347;
t327 = t271 * t348;
t250 = pkin(2) * t327;
t166 = t250 + t283;
t68 = -qJD(4) * t116 - t166 * t274 + t199 * t277;
t43 = pkin(4) * t326 - qJ(5) * t184 + qJD(5) * t291 + t68;
t185 = qJD(4) * t291 + t271 * t324;
t67 = t166 * t277 + t172 * t344 - t190 * t345 + t199 * t274;
t47 = qJ(5) * t185 + qJD(5) * t228 + t67;
t17 = t270 * t43 + t361 * t47;
t115 = t172 * t277 - t190 * t274;
t89 = pkin(4) * t358 + qJ(5) * t291 + t115;
t96 = qJ(5) * t228 + t116;
t49 = t270 * t89 + t361 * t96;
t393 = mrSges(5,3) * t207;
t392 = mrSges(5,3) * t208;
t383 = t102 * Ifges(6,4);
t382 = t103 * Ifges(6,1);
t381 = t103 * Ifges(6,4);
t137 = -mrSges(5,1) * t207 + mrSges(5,2) * t208;
t218 = -mrSges(4,1) * t328 - mrSges(4,3) * t257;
t354 = -t218 + t137;
t353 = t257 * t467 + t329 * t468;
t231 = pkin(8) * t357 + t261;
t342 = qJD(6) * t273;
t341 = qJD(6) * t276;
t339 = t272 * t398;
t335 = -0.3e1 / 0.2e1 * Ifges(4,6) - 0.3e1 / 0.2e1 * Ifges(3,4);
t211 = -qJ(3) * t272 - t231;
t330 = t361 * pkin(4);
t63 = mrSges(6,1) * t102 + mrSges(6,2) * t103;
t319 = t427 * t358;
t189 = pkin(3) * t357 - t211;
t243 = qJD(2) * t253;
t245 = t257 * qJD(3);
t295 = qJD(2) * t319;
t149 = -qJD(1) * t295 + t243 + t245;
t104 = -pkin(4) * t158 + t149;
t29 = pkin(5) * t102 - pkin(10) * t103 + t104;
t6 = pkin(10) * t316 + t8;
t1 = qJD(6) * t9 + t273 * t29 + t276 * t6;
t2 = -qJD(6) * t10 - t273 * t6 + t276 * t29;
t315 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t314 = t1 * t273 + t2 * t276;
t313 = t10 * t276 - t273 * t9;
t311 = -mrSges(7,1) * t276 + mrSges(7,2) * t273;
t307 = Ifges(7,1) * t273 + t385;
t304 = Ifges(7,2) * t276 + t386;
t301 = Ifges(7,5) * t273 + Ifges(7,6) * t276;
t45 = pkin(10) * t358 + t49;
t141 = -pkin(4) * t228 + t189;
t72 = pkin(5) * t150 - pkin(10) * t151 + t141;
t22 = t273 * t72 + t276 * t45;
t21 = -t273 * t45 + t276 * t72;
t73 = -mrSges(7,2) * t128 + mrSges(7,3) * t111;
t74 = mrSges(7,1) * t128 - mrSges(7,3) * t112;
t298 = -t273 * t74 + t276 * t73;
t134 = mrSges(5,1) * t316 - mrSges(5,3) * t157;
t135 = -mrSges(5,2) * t316 + mrSges(5,3) * t158;
t294 = t277 * t134 + t274 * t135;
t162 = -mrSges(5,2) * t244 + t393;
t163 = mrSges(5,1) * t244 - t392;
t293 = t277 * t162 - t274 * t163;
t254 = qJD(2) * t339;
t224 = -pkin(8) * t327 + t254;
t117 = -mrSges(6,2) * t244 + mrSges(6,3) * t320;
t292 = -t117 - t298;
t126 = -t151 * t273 + t276 * t358;
t127 = t151 * t276 + t273 * t358;
t16 = -t270 * t47 + t361 * t43;
t48 = -t270 * t96 + t361 * t89;
t209 = -pkin(8) * t317 + t243;
t267 = t272 * qJD(3);
t171 = t254 + t267 - t295;
t288 = (-qJ(3) * t347 - t346) * t271;
t225 = t231 * qJD(2);
t179 = -t209 - t245;
t210 = qJD(1) * t225;
t286 = -t209 * mrSges(3,2) - t179 * mrSges(4,3) + t210 * t467;
t119 = -pkin(4) * t185 + t171;
t282 = t1 * t276 - t2 * t273 + (-t10 * t273 - t276 * t9) * qJD(6);
t264 = -t330 - pkin(5);
t241 = Ifges(3,5) * t316;
t240 = Ifges(4,5) * t317;
t239 = Ifges(5,3) * t316;
t238 = Ifges(6,3) * t316;
t230 = -t258 + t339;
t221 = -qJ(3) * t328 + t248;
t220 = (mrSges(4,2) * t278 - mrSges(4,3) * t275) * t351;
t217 = -mrSges(3,2) * t257 + mrSges(3,3) * t328;
t213 = t272 * t331 + t258;
t205 = -t224 - t267;
t200 = t250 + t288;
t197 = -qJD(1) * t319 + t253;
t176 = -t236 * t270 + t321 * t361;
t175 = qJD(1) * t288 + t242;
t145 = -t187 * t276 + t257 * t273;
t144 = t187 * t273 + t257 * t276;
t114 = t184 * t361 + t185 * t270;
t113 = t184 * t270 - t185 * t361;
t105 = -mrSges(5,1) * t158 + mrSges(5,2) * t157;
t92 = Ifges(5,1) * t157 + Ifges(5,4) * t158 + Ifges(5,5) * t316;
t91 = Ifges(5,4) * t157 + Ifges(5,2) * t158 + Ifges(5,6) * t316;
t87 = -mrSges(6,2) * t316 - mrSges(6,3) * t102;
t84 = -mrSges(6,1) * t320 + mrSges(6,2) * t289;
t71 = pkin(5) * t289 - pkin(10) * t320 + t397;
t70 = -qJD(6) * t127 - t114 * t273 + t276 * t326;
t69 = qJD(6) * t126 + t114 * t276 + t273 * t326;
t51 = Ifges(6,5) * t316 + t382 - t383;
t50 = -t102 * Ifges(6,2) + Ifges(6,6) * t316 + t381;
t44 = -pkin(5) * t358 - t48;
t38 = pkin(5) * t113 - pkin(10) * t114 + t119;
t37 = t361 * t82 - t366;
t36 = t270 * t82 + t80;
t19 = t273 * t71 + t276 * t37;
t18 = -t273 * t37 + t276 * t71;
t15 = pkin(10) * t326 + t17;
t14 = -pkin(5) * t326 - t16;
t13 = Ifges(7,1) * t58 + Ifges(7,4) * t59 + Ifges(7,5) * t102;
t12 = Ifges(7,4) * t58 + Ifges(7,2) * t59 + Ifges(7,6) * t102;
t5 = -pkin(5) * t316 - t7;
t4 = -qJD(6) * t22 - t15 * t273 + t276 * t38;
t3 = qJD(6) * t21 + t15 * t276 + t273 * t38;
t23 = [((-mrSges(4,1) * t179 + mrSges(4,2) * t175 + mrSges(3,3) * t209) * t278 + (-t175 * mrSges(4,3) + t238 / 0.2e1 + t239 / 0.2e1 + t468 * t210 + t442) * t275) * t271 + (t275 * t444 + t278 * t441) * t349 + m(3) * (t209 * t231 - t210 * t230 + t222 * t225 + t223 * t224) + m(4) * (t175 * t212 + t179 * t211 + t186 * t225 + t195 * t205 + t200 * t203 + t210 * t213) + m(5) * (t115 * t61 + t116 * t60 + t149 * t189 + t161 * t171 + t67 * t94 + t68 * t93) + m(7) * (t1 * t22 + t10 * t3 + t14 * t27 + t2 * t21 + t4 * t9 + t44 * t5) + m(6) * (t104 * t141 + t119 * t122 + t16 * t32 + t17 * t33 + t48 * t7 + t49 * t8) + t207 * (Ifges(5,4) * t184 + Ifges(5,2) * t185) / 0.2e1 + t113 * t435 + t150 * t437 + t92 * t406 + t91 * t407 + (Ifges(5,1) * t184 + Ifges(5,4) * t185) * t408 + t51 * t413 + t50 * t414 + (Ifges(7,5) * t69 + Ifges(7,6) * t70 + Ifges(7,3) * t113) * t417 + t320 * (Ifges(6,4) * t114 - Ifges(6,2) * t113) / 0.2e1 + t289 * (Ifges(6,1) * t114 - Ifges(6,4) * t113) / 0.2e1 + (-t113 * t33 - t114 * t32 - t150 * t8 - t151 * t7) * mrSges(6,3) + t353 * t225 + (Ifges(7,1) * t69 + Ifges(7,4) * t70 + Ifges(7,5) * t113) * t419 + (Ifges(7,4) * t69 + Ifges(7,2) * t70 + Ifges(7,6) * t113) * t421 + (Ifges(7,5) * t127 + Ifges(7,6) * t126 + Ifges(7,3) * t150) * t423 + (Ifges(6,4) * t151 - Ifges(6,2) * t150) * t424 + t113 * t430 + (Ifges(7,4) * t127 + Ifges(7,2) * t126 + Ifges(7,6) * t150) * t431 + (Ifges(7,1) * t127 + Ifges(7,4) * t126 + Ifges(7,5) * t150) * t432 + t69 * t433 + (-Ifges(5,4) * t291 + Ifges(5,2) * t228) * t411 + (-Ifges(5,1) * t291 + Ifges(5,4) * t228) * t412 + (-t184 * t93 + t185 * t94 + t228 * t60 + t291 * t61) * mrSges(5,3) + t149 * (-mrSges(5,1) * t228 - mrSges(5,2) * t291) + ((t213 * mrSges(4,1) - t230 * mrSges(3,3) - t212 * mrSges(4,3) + Ifges(5,5) * t406 + Ifges(5,6) * t407 + Ifges(6,5) * t413 + Ifges(6,6) * t414 + (-Ifges(4,4) + t438) * t272 + (-t278 * t335 - 0.2e1 * t425) * t271) * t278 + (t211 * mrSges(4,1) - t212 * mrSges(4,2) - t231 * mrSges(3,3) + (-Ifges(3,6) + Ifges(4,5) / 0.2e1) * t272 + (t275 * t335 - 0.2e1 * t426) * t271 + (-0.3e1 / 0.2e1 * Ifges(4,3) + 0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2) + 0.3e1 / 0.2e1 * Ifges(4,2) + t334) * t357) * t275) * t323 + t224 * t217 + t200 * t220 + (t240 / 0.2e1 + t241 / 0.2e1 + t286) * t272 + t205 * t218 + t189 * t105 + t184 * t125 / 0.2e1 + t161 * (-mrSges(5,1) * t185 + mrSges(5,2) * t184) + t185 * t124 / 0.2e1 + t67 * t162 + t68 * t163 + t171 * t137 + t1 * (-mrSges(7,2) * t150 + mrSges(7,3) * t126) + t2 * (mrSges(7,1) * t150 - mrSges(7,3) * t127) + t104 * (mrSges(6,1) * t150 + mrSges(6,2) * t151) + t115 * t134 + t116 * t135 + t141 * t63 + t5 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t127 * t13 / 0.2e1 + t119 * t84 + t122 * (mrSges(6,1) * t113 + mrSges(6,2) * t114) + t126 * t12 / 0.2e1 + t114 * t79 / 0.2e1 + t17 * t117 + t16 * t118 + t9 * (mrSges(7,1) * t113 - mrSges(7,3) * t69) + t10 * (-mrSges(7,2) * t113 + mrSges(7,3) * t70) + (Ifges(5,5) * t184 + Ifges(6,5) * t114 + Ifges(5,6) * t185 - Ifges(6,6) * t113) * t453 + t49 * t87 + t48 * t88 + t3 * t73 + t4 * t74 + t27 * (-mrSges(7,1) * t70 + mrSges(7,2) * t69) + t70 * t41 / 0.2e1 + t14 * t66 + t44 * t20 + t22 * t26 + t21 * t25 + t103 * (Ifges(6,1) * t151 - Ifges(6,4) * t150) / 0.2e1; (t430 + t435 + t454 + t469) * t187 + t465 * t74 + (t1 * t108 + t10 * t466 + t107 * t2 + t176 * t5 + t27 * t463 + t465 * t9) * m(7) + t466 * t73 + (-t355 * mrSges(7,2) + t227 * t391) * t10 + t462 * t117 + t463 * t66 + (t104 * t265 + t122 * t461 - t176 * t7 + t177 * t8 + t32 * t464 + t33 * t462) * m(6) + t464 * t118 + (-pkin(2) * t210 - qJ(3) * t179 - t186 * t223 + t195 * t458 - t203 * t221) * m(4) + t459 * mrSges(6,3) + t461 * t84 + (t27 * mrSges(7,1) + Ifges(7,4) * t420 + Ifges(7,2) * t422 + Ifges(7,6) * t418 + t456) * (-t188 * t273 + t276 * t328) + (-t368 / 0.2e1 - t377 / 0.2e1 - t374 / 0.2e1 + t457) * t227 + (t428 + t455) * t188 + (-t378 / 0.2e1 - t379 / 0.2e1 - t380 / 0.2e1 + t367 / 0.2e1 + t376 / 0.2e1 + t373 / 0.2e1 + t443) * t226 + ((-t247 / 0.2e1 + (t425 - t384 / 0.2e1) * t351 - t441) * t278 + (((t426 + (Ifges(4,6) / 0.2e1 + Ifges(3,4) / 0.2e1) * t275) * t271 + (-Ifges(4,2) / 0.2e1 + Ifges(4,3) / 0.2e1 - Ifges(3,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t357) * qJD(1) + t439 - t444) * t275 + (t300 * mrSges(4,1) - Ifges(4,4) * t278 - Ifges(3,6) * t275 + (Ifges(5,5) * t277 - Ifges(6,5) * t233 - Ifges(5,6) * t274 - Ifges(6,6) * t234) * t278 / 0.2e1) * qJD(2)) * t351 + ((-m(5) * t296 + t293) * t279 + t439) * qJD(4) + t240 + t241 + (-t5 * t310 + t302 * t424 - t59 * t305 / 0.2e1 - t58 * t308 / 0.2e1 - t382 / 0.2e1 - t104 * mrSges(6,2) + t383 / 0.2e1 - t51 / 0.2e1 + t13 * t401 + t12 * t403 + t7 * mrSges(6,3) + t314 * mrSges(7,3) + (mrSges(7,3) * t313 + t27 * t311 + t301 * t417 + t304 * t421 + t307 * t419 + t400 * t41 + t403 * t42) * qJD(6)) * t233 + m(5) * (qJ(3) * t149 + qJD(3) * t161 + t279 * t446) + (-Ifges(5,2) * t274 + t388) * t411 + (Ifges(5,1) * t277 - t389) * t412 - m(5) * (t120 * t93 + t121 * t94 + t161 * t197) - t353 * t223 + t354 * qJD(3) + (mrSges(6,1) * t355 - mrSges(6,2) * t352) * t122 + (Ifges(7,5) * t418 + Ifges(7,1) * t420 + Ifges(7,4) * t422 + t9 * mrSges(7,3) - t27 * mrSges(7,2) - t42 / 0.2e1) * (t188 * t276 + t273 * t328) + t286 + (-t218 + t217) * t222 - t446 * mrSges(5,3) + t394 * t176 + (-t8 * mrSges(6,3) - t381 / 0.2e1 + t57 / 0.2e1 + t56 / 0.2e1 + t437 - t50 / 0.2e1 + t104 * mrSges(6,1) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t102 + t315) * t234 + t149 * (mrSges(5,1) * t274 + mrSges(5,2) * t277) + t277 * t92 / 0.2e1 + t265 * t63 - t221 * t220 - t197 * t137 + t177 * t87 - t121 * t162 - t120 * t163 + qJ(3) * t105 + t107 * t25 + t108 * t26 + t294 * t279 + t91 * t402; t187 * t117 - t144 * t74 - t145 * t73 + t394 * t233 + t293 * qJD(4) + (-t84 - t354) * t257 + t292 * t226 + (mrSges(4,1) * t347 + (t220 + t293) * t275) * t351 + (t87 + (-t273 * t73 - t276 * t74) * qJD(6) + t445) * t234 + t294 + t447 * t352 + (-t10 * t145 - t144 * t9 - t226 * t313 + t233 * t5 + t234 * t282 + t27 * t352) * m(7) + (-t122 * t257 - t233 * t7 + t234 * t8 - t459) * m(6) + (-t161 * t257 - t244 * t296 + t446) * m(5) + (t195 * t257 + t203 * t329 + t210) * m(4); t460 * qJD(6) + t456 * t342 + (-t122 * mrSges(6,2) + t32 * mrSges(6,3) + t10 * t391 + t455 + t457) * t320 + (-mrSges(6,1) * t122 + mrSges(7,2) * t10 + mrSges(6,3) * t33 + t443 - t454) * t289 + t442 + t238 + t239 + (t163 + t392) * t94 + t5 * t311 + (m(7) * t282 - t341 * t74 - t342 * t73 + t445) * (pkin(10) + t396) + t13 * t403 + (Ifges(5,5) * t207 - Ifges(5,6) * t208) * t405 + t124 * t408 + (Ifges(5,1) * t207 - t370) * t409 + (-t162 + t393) * t93 + t301 * t423 + t304 * t431 + t307 * t432 + t341 * t433 + (-t122 * t397 + t32 * t36 - t33 * t37 + (t270 * t8 + t361 * t7) * pkin(4)) * m(6) - t84 * t397 - t447 * t36 + (t111 * t305 + t112 * t308 + t128 * t302) * qJD(6) / 0.2e1 + t88 * t330 - t2 * t391 + (-Ifges(5,2) * t208 + t125 + t204) * t410 + t264 * t20 - t161 * (mrSges(5,1) * t208 + mrSges(5,2) * t207) + (-t10 * t19 - t18 * t9 + t264 * t5 - t27 * t36) * m(7) - t37 * t117 - t19 * t73 - t18 * t74 + t1 * t390 + t87 * t396 + t12 * t400; t276 * t25 + t273 * t26 - t447 * t289 + t298 * qJD(6) + t292 * t320 + t63 + (t128 * t313 - t27 * t289 + t314) * m(7) + (t289 * t32 - t320 * t33 + t104) * m(6); -t27 * (mrSges(7,1) * t112 + mrSges(7,2) * t111) + (Ifges(7,1) * t111 - t387) * t420 + t41 * t419 + (Ifges(7,5) * t111 - Ifges(7,6) * t112) * t418 - t9 * t73 + t10 * t74 + (t10 * t112 + t111 * t9) * mrSges(7,3) + t315 + t11 + (-Ifges(7,2) * t112 + t110 + t42) * t422;];
tauc  = t23(:);
