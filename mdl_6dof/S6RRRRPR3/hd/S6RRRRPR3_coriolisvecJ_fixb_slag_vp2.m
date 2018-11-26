% Calculate vector of centrifugal and coriolis load on the joints for
% S6RRRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
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
% Datum: 2018-11-23 18:13
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RRRRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 18:13:05
% EndTime: 2018-11-23 18:13:16
% DurationCPUTime: 11.59s
% Computational Cost: add. (14938->606), mult. (38216->807), div. (0->0), fcn. (27953->8), ass. (0->302)
t462 = Ifges(6,3) + Ifges(5,2);
t392 = sin(qJ(3));
t393 = sin(qJ(2));
t395 = cos(qJ(3));
t396 = cos(qJ(2));
t294 = t392 * t393 - t395 * t396;
t394 = cos(qJ(4));
t286 = t394 * t294;
t204 = qJD(1) * t286;
t231 = t392 * t396 + t393 * t395;
t220 = t231 * qJD(1);
t260 = sin(qJ(4));
t175 = t220 * t260 + t204;
t405 = t175 / 0.2e1;
t258 = qJD(2) + qJD(3);
t257 = qJD(4) + t258;
t399 = -t257 / 0.2e1;
t398 = t257 / 0.2e1;
t290 = qJD(1) * t294;
t329 = qJD(1) * t396;
t195 = -qJD(1) * pkin(1) - pkin(2) * t329 + pkin(3) * t290;
t285 = t260 * t290;
t282 = t220 * t394 - t285;
t100 = t175 * pkin(4) - qJ(5) * t282 + t195;
t440 = Ifges(6,6) * t282;
t441 = Ifges(5,4) * t282;
t341 = t396 * pkin(7);
t244 = pkin(8) * t396 + t341;
t236 = t244 * qJD(1);
t224 = t395 * t236;
t338 = t393 * pkin(7);
t243 = -pkin(8) * t393 - t338;
t235 = t243 * qJD(1);
t371 = qJD(2) * pkin(2);
t227 = t235 + t371;
t187 = t227 * t392 + t224;
t288 = pkin(9) * t290;
t150 = -t288 + t187;
t142 = t394 * t150;
t221 = t392 * t236;
t186 = t395 * t227 - t221;
t214 = t220 * pkin(9);
t149 = -t214 + t186;
t295 = t258 * pkin(3) + t149;
t93 = t260 * t295 + t142;
t77 = -t257 * qJ(5) - t93;
t461 = t195 * mrSges(5,1) + t77 * mrSges(6,1) - t100 * mrSges(6,2) - t93 * mrSges(5,3) + Ifges(6,5) * t398 - t440 / 0.2e1 + Ifges(5,6) * t399 - t441 / 0.2e1 + t462 * t405;
t460 = mrSges(6,2) - mrSges(5,1);
t459 = mrSges(5,3) + mrSges(6,1);
t324 = qJD(4) * t394;
t352 = t260 * t150;
t99 = t149 * t394 - t352;
t458 = pkin(3) * t324 - t99;
t358 = qJ(5) * t175;
t126 = pkin(4) * t282 + t358;
t259 = sin(qJ(6));
t261 = cos(qJ(6));
t307 = Ifges(7,5) * t259 + Ifges(7,6) * t261;
t373 = Ifges(7,4) * t259;
t308 = Ifges(7,2) * t261 + t373;
t372 = Ifges(7,4) * t261;
t309 = Ifges(7,1) * t259 + t372;
t397 = -t259 / 0.2e1;
t170 = qJD(6) + t282;
t407 = -t170 / 0.2e1;
t152 = t175 * t259 + t257 * t261;
t409 = -t152 / 0.2e1;
t151 = t175 * t261 - t257 * t259;
t410 = -t151 / 0.2e1;
t310 = mrSges(7,1) * t261 - mrSges(7,2) * t259;
t455 = t175 * pkin(5);
t54 = -t77 - t455;
t374 = Ifges(7,4) * t152;
t63 = t151 * Ifges(7,2) + t170 * Ifges(7,6) + t374;
t449 = t54 * t310 - t261 * t63 / 0.2e1;
t146 = Ifges(7,4) * t151;
t64 = t152 * Ifges(7,1) + t170 * Ifges(7,5) + t146;
t457 = t307 * t407 + t308 * t410 + t309 * t409 + t64 * t397 + t449 - t461;
t454 = -Ifges(5,5) + Ifges(6,4);
t453 = -Ifges(5,6) + Ifges(6,5);
t166 = Ifges(5,4) * t175;
t120 = Ifges(5,1) * t282 + t257 * Ifges(5,5) - t166;
t365 = t170 * Ifges(7,3);
t366 = t152 * Ifges(7,5);
t367 = t151 * Ifges(7,6);
t62 = t365 + t366 + t367;
t452 = t120 + t62;
t436 = -qJD(5) - t458;
t193 = -t235 * t392 - t224;
t153 = t288 + t193;
t194 = t395 * t235 - t221;
t154 = -t214 + t194;
t103 = t260 * t153 + t154 * t394;
t340 = t395 * pkin(2);
t252 = t340 + pkin(3);
t325 = qJD(3) * t395;
t319 = pkin(2) * t325;
t323 = t392 * qJD(3);
t182 = -(qJD(4) * t392 + t323) * pkin(2) * t260 + t252 * t324 + t394 * t319;
t181 = -qJD(5) - t182;
t451 = t181 + t103;
t275 = t258 * t231;
t269 = t394 * t275;
t274 = t258 * t294;
t272 = qJD(1) * t274;
t91 = qJD(1) * t269 - qJD(4) * t285 + t220 * t324 - t260 * t272;
t48 = qJD(6) * t151 + t259 * t91;
t268 = t394 * t274;
t273 = qJD(1) * t275;
t347 = qJD(4) * t260;
t90 = qJD(1) * t268 + qJD(4) * t204 + t220 * t347 + t260 * t273;
t26 = -mrSges(7,1) * t90 - mrSges(7,3) * t48;
t49 = -qJD(6) * t152 + t261 * t91;
t27 = mrSges(7,2) * t90 + mrSges(7,3) * t49;
t237 = t243 * qJD(2);
t228 = qJD(1) * t237;
t238 = t244 * qJD(2);
t229 = qJD(1) * t238;
t299 = -t228 * t392 - t229 * t395;
t266 = pkin(9) * t272 - t227 * t323 - t236 * t325 + t299;
t132 = t227 * t325 + t395 * t228 - t392 * t229 - t236 * t323;
t97 = -pkin(9) * t273 + t132;
t21 = qJD(4) * t93 + t260 * t97 - t394 * t266;
t11 = -t90 * pkin(5) + t21;
t326 = qJD(2) * t393;
t315 = qJD(1) * t326;
t247 = pkin(2) * t315;
t159 = pkin(3) * t273 + t247;
t23 = t91 * pkin(4) + t90 * qJ(5) - qJD(5) * t282 + t159;
t12 = t91 * pkin(10) + t23;
t415 = pkin(4) + pkin(10);
t387 = t282 * pkin(5);
t137 = t394 * t295;
t92 = -t137 + t352;
t303 = t92 + t387;
t448 = qJD(5) + t303;
t51 = -t257 * t415 + t448;
t65 = t175 * pkin(10) + t100;
t24 = -t259 * t65 + t261 * t51;
t1 = qJD(6) * t24 + t11 * t259 + t12 * t261;
t25 = t259 * t51 + t261 * t65;
t357 = qJD(6) * t25;
t2 = t11 * t261 - t12 * t259 - t357;
t305 = t24 * t259 - t25 * t261;
t283 = m(7) * (-qJD(6) * t305 + t1 * t259 + t2 * t261);
t450 = t259 * t27 + t261 * t26 + t283;
t434 = -qJD(5) - t92;
t76 = -pkin(4) * t257 - t434;
t446 = t76 * mrSges(6,1) + t24 * mrSges(7,1) + t195 * mrSges(5,2) - t25 * mrSges(7,2) + t92 * mrSges(5,3) - t100 * mrSges(6,3);
t447 = -t446 - t365 / 0.2e1 - t367 / 0.2e1;
t418 = t48 / 0.2e1;
t417 = t49 / 0.2e1;
t416 = -t90 / 0.2e1;
t444 = m(5) * t92;
t16 = -mrSges(7,1) * t49 + mrSges(7,2) * t48;
t438 = -t91 * mrSges(6,1) + t16;
t437 = t387 - t451;
t435 = t387 - t436;
t155 = mrSges(6,1) * t175 - mrSges(6,3) * t257;
t157 = -mrSges(5,2) * t257 - mrSges(5,3) * t175;
t433 = t155 - t157;
t432 = -t460 * t257 - t459 * t282;
t328 = qJD(1) * t393;
t431 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t328) * t341 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t329) * t338;
t427 = m(6) * t76 - t432;
t426 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t292 = Ifges(4,4) * t294;
t425 = -(-Ifges(4,2) * t231 - t292) * t290 + t258 * (-Ifges(4,5) * t294 - Ifges(4,6) * t231) / 0.2e1;
t108 = -mrSges(7,2) * t170 + mrSges(7,3) * t151;
t109 = mrSges(7,1) * t170 - mrSges(7,3) * t152;
t345 = qJD(6) * t261;
t346 = qJD(6) * t259;
t424 = t108 * t345 - t109 * t346 + t450;
t422 = m(7) * (t24 * t261 + t25 * t259) + t261 * t109 + t259 * t108 + t427;
t419 = Ifges(7,1) * t418 + Ifges(7,4) * t417 + Ifges(7,5) * t416;
t165 = Ifges(6,6) * t175;
t118 = t257 * Ifges(6,4) - Ifges(6,2) * t282 + t165;
t412 = t118 / 0.2e1;
t408 = t152 / 0.2e1;
t406 = -t175 / 0.2e1;
t403 = -t282 / 0.2e1;
t402 = t282 / 0.2e1;
t400 = t220 / 0.2e1;
t390 = pkin(3) * t220;
t389 = pkin(3) * t260;
t383 = t90 * mrSges(6,1);
t382 = t90 * mrSges(5,3);
t380 = t91 * mrSges(5,3);
t376 = Ifges(4,4) * t220;
t375 = Ifges(4,4) * t231;
t196 = t395 * t243 - t244 * t392;
t167 = -t231 * pkin(9) + t196;
t197 = t392 * t243 + t395 * t244;
t168 = -pkin(9) * t294 + t197;
t121 = -t394 * t167 + t168 * t260;
t370 = t121 * t21;
t363 = t220 * mrSges(4,3);
t362 = t259 * mrSges(7,3);
t360 = t261 * mrSges(7,3);
t291 = t260 * t294;
t106 = -qJD(4) * t291 + t231 * t324 - t260 * t274 + t269;
t356 = t106 * t259;
t191 = t231 * t260 + t286;
t355 = t191 * t259;
t354 = t191 * t261;
t351 = t261 * t106;
t101 = -mrSges(7,1) * t151 + mrSges(7,2) * t152;
t349 = -t155 + t101;
t253 = -pkin(2) * t396 - pkin(1);
t242 = qJD(1) * t253;
t344 = qJD(1) * qJD(2);
t343 = Ifges(7,5) * t48 + Ifges(7,6) * t49 - Ifges(7,3) * t90;
t342 = pkin(2) * t392;
t339 = t394 * pkin(3);
t336 = Ifges(3,4) * t396;
t335 = Ifges(3,4) * t393;
t327 = qJD(2) * t396;
t322 = -t346 / 0.2e1;
t98 = t149 * t260 + t142;
t102 = -t394 * t153 + t154 * t260;
t256 = pkin(2) * t326;
t255 = pkin(2) * t328;
t251 = -t339 - pkin(4);
t317 = t394 * t392;
t316 = qJD(1) * t327;
t313 = t327 / 0.2e1;
t311 = -t237 * t392 - t395 * t238;
t215 = t252 * t394 - t260 * t342;
t192 = t231 * t394 - t291;
t202 = pkin(3) * t294 + t253;
t281 = -t192 * qJ(5) + t202;
t75 = t191 * t415 + t281;
t78 = pkin(5) * t192 + t121;
t38 = t259 * t78 + t261 * t75;
t37 = -t259 * t75 + t261 * t78;
t304 = t108 * t261 - t109 * t259;
t211 = -pkin(4) - t215;
t122 = t260 * t167 + t168 * t394;
t302 = t191 * t345 + t356;
t301 = t191 * t346 - t351;
t20 = qJD(4) * t137 - t150 * t347 + t260 * t266 + t394 * t97;
t138 = t395 * t237 - t392 * t238 + t243 * t325 - t244 * t323;
t113 = -pkin(9) * t275 + t138;
t267 = pkin(9) * t274 - t243 * t323 - t244 * t325 + t311;
t35 = -t394 * t113 - t167 * t324 + t168 * t347 - t260 * t267;
t300 = Ifges(3,5) * t396 - Ifges(3,6) * t393;
t107 = t390 + t126;
t216 = pkin(2) * t317 + t260 * t252;
t104 = t107 + t255;
t297 = pkin(1) * (mrSges(3,1) * t393 + mrSges(3,2) * t396);
t296 = (Ifges(3,2) * t396 + t335) * qJD(1);
t17 = -qJD(5) * t257 - t20;
t293 = (Ifges(3,1) * t396 - t335) * t393;
t36 = qJD(4) * t122 + t260 * t113 - t394 * t267;
t287 = mrSges(4,3) * t290;
t10 = -pkin(5) * t91 - t17;
t8 = t48 * Ifges(7,4) + t49 * Ifges(7,2) - t90 * Ifges(7,6);
t280 = -t20 * mrSges(5,2) - t17 * mrSges(6,3) + t24 * mrSges(7,3) * t346 + (Ifges(7,1) * t261 - t373) * t418 + (-Ifges(7,2) * t259 + t372) * t417 + t10 * (mrSges(7,1) * t259 + mrSges(7,2) * t261) + t64 * t322 + t8 * t397 + t261 * t419 + (Ifges(7,5) * t261 - Ifges(7,6) * t259) * t416 + t453 * t91 + t454 * t90 + t460 * t21 + t449 * qJD(6) - (t151 * t308 + t152 * t309 + t170 * t307) * qJD(6) / 0.2e1;
t271 = mrSges(4,3) * t273;
t270 = mrSges(4,3) * t272;
t178 = pkin(3) * t275 + t256;
t105 = qJD(4) * t286 + t231 * t347 + t260 * t275 + t268;
t34 = t106 * pkin(4) + t105 * qJ(5) - t192 * qJD(5) + t178;
t265 = t258 * (-Ifges(4,1) * t294 - t375);
t133 = -qJD(3) * t187 + t299;
t163 = -Ifges(4,2) * t290 + Ifges(4,6) * t258 + t376;
t213 = Ifges(4,4) * t290;
t164 = Ifges(4,1) * t220 + Ifges(4,5) * t258 - t213;
t263 = -t242 * (t220 * mrSges(4,1) - mrSges(4,2) * t290) - t258 * (-Ifges(4,5) * t290 - Ifges(4,6) * t220) / 0.2e1 - t220 * (-Ifges(4,1) * t290 - t376) / 0.2e1 - t25 * mrSges(7,3) * t345 - t441 * t403 + t440 * t402 - Ifges(4,5) * t272 - Ifges(4,6) * t273 + t280 - t186 * t287 - t132 * mrSges(4,2) + t133 * mrSges(4,1) - t1 * t362 - t2 * t360 + t187 * t363 + t163 * t400 + t452 * t405 + (-Ifges(4,2) * t220 + t164 - t213) * t290 / 0.2e1 + (-Ifges(5,1) * t403 - Ifges(5,4) * t405 - Ifges(7,5) * t409 + Ifges(6,2) * t402 + Ifges(6,6) * t406 - Ifges(7,6) * t410 - Ifges(7,3) * t407 + t399 * t454 - t412 + t446) * t175 + (-Ifges(5,2) * t405 + Ifges(6,3) * t406 + t24 * t362 - t360 * t25 + t399 * t453 + t457) * t282;
t254 = Ifges(3,4) * t329;
t250 = qJ(5) + t389;
t219 = Ifges(3,1) * t328 + Ifges(3,5) * qJD(2) + t254;
t218 = Ifges(3,6) * qJD(2) + t296;
t210 = qJ(5) + t216;
t200 = mrSges(4,1) * t258 - t363;
t199 = -t258 * mrSges(4,2) - t287;
t198 = t255 + t390;
t185 = mrSges(4,1) * t290 + t220 * mrSges(4,2);
t169 = t282 * pkin(10);
t139 = -qJD(3) * t197 + t311;
t130 = -mrSges(6,2) * t175 - mrSges(6,3) * t282;
t129 = mrSges(5,1) * t175 + mrSges(5,2) * t282;
t114 = t191 * pkin(4) + t281;
t89 = t282 * t415 + t358;
t79 = -t191 * pkin(5) + t122;
t69 = t107 + t169;
t68 = t104 + t169;
t66 = t102 - t455;
t60 = t98 - t455;
t57 = t93 - t455;
t33 = t259 * t57 + t261 * t89;
t32 = -t259 * t89 + t261 * t57;
t31 = t259 * t66 + t261 * t68;
t30 = -t259 * t68 + t261 * t66;
t29 = t259 * t60 + t261 * t69;
t28 = -t259 * t69 + t261 * t60;
t15 = -t105 * pkin(5) + t36;
t14 = -pkin(5) * t106 - t35;
t13 = t106 * pkin(10) + t34;
t4 = -qJD(6) * t38 - t13 * t259 + t15 * t261;
t3 = qJD(6) * t37 + t13 * t261 + t15 * t259;
t5 = [-(Ifges(4,1) * t231 - t292) * t272 / 0.2e1 - (t296 + t218) * t326 / 0.2e1 + (qJD(6) * t64 + t8) * t354 / 0.2e1 + (-0.2e1 * t297 + t293) * t344 + t459 * (-t121 * t90 - t122 * t91 + t192 * t21) + (-Ifges(5,4) * t402 - Ifges(5,2) * t406 + Ifges(6,3) * t405 + t453 * t398 + t461) * t106 + ((Ifges(3,1) * t393 + t336) * t313 + t231 * t265 / 0.2e1) * qJD(1) + t170 * (Ifges(7,5) * t302 - Ifges(7,6) * t301) / 0.2e1 + t151 * (Ifges(7,4) * t302 - Ifges(7,2) * t301) / 0.2e1 + (Ifges(7,1) * t302 - Ifges(7,4) * t301) * t408 + (Ifges(6,6) * t405 - t452 / 0.2e1 - Ifges(5,1) * t402 - Ifges(5,4) * t406 - Ifges(7,5) * t408 + Ifges(6,2) * t403 + t454 * t398 + t412 + t447) * t105 + (t1 * t354 - t2 * t355 - t24 * t302 - t25 * t301) * mrSges(7,3) + t425 * qJD(3) + 0.2e1 * t242 * t258 * (mrSges(4,1) * t231 - mrSges(4,2) * t294) + (t106 * t403 + t191 * t90) * Ifges(6,6) + m(5) * (t122 * t20 + t159 * t202 + t178 * t195 + t370) + m(6) * (t100 * t34 + t114 * t23 - t122 * t17 + t370) + (t159 * mrSges(5,1) + t17 * mrSges(6,1) - t23 * mrSges(6,2) - t20 * mrSges(5,3) - t10 * t310 + t307 * t416 + t308 * t417 + t309 * t418 + t63 * t322 + t462 * t91 + (t90 / 0.2e1 - t416) * Ifges(5,4)) * t191 + (t427 + t444) * t36 + (-m(5) * t93 + m(6) * t77 + t433) * t35 + (t300 * qJD(2) / 0.2e1 + t425 - t431) * qJD(2) + t185 * t256 + t219 * t313 + t196 * t270 - t275 * t163 / 0.2e1 - t274 * t164 / 0.2e1 + m(7) * (t1 * t38 + t10 * t79 + t14 * t54 + t2 * t37 + t24 * t4 + t25 * t3) + t54 * (mrSges(7,1) * t301 + mrSges(7,2) * t302) + (mrSges(4,1) * t294 + t231 * mrSges(4,2)) * t247 - t197 * t271 + (-t132 * t294 - t133 * t231 + t186 * t274 - t187 * t275) * mrSges(4,3) + t138 * t199 + t139 * t200 + t202 * (mrSges(5,1) * t91 - mrSges(5,2) * t90) + t178 * t129 + t34 * t130 + m(4) * (t132 * t197 + t133 * t196 + t187 * t138 + t186 * t139 + (t242 * t393 + t253 * t328) * t371) + (-Ifges(3,2) * t393 + t336) * t316 + t114 * (-mrSges(6,2) * t91 + mrSges(6,3) * t90) + t3 * t108 + t4 * t109 + t14 * t101 - (-Ifges(4,2) * t294 + t375) * t273 / 0.2e1 + t79 * t16 + t64 * t356 / 0.2e1 + t63 * t351 / 0.2e1 + t38 * t27 + t37 * t26 + (t159 * mrSges(5,2) - t23 * mrSges(6,3) + Ifges(7,6) * t417 + Ifges(7,5) * t418 + (Ifges(7,3) + Ifges(5,1)) * t416 + t426 + t343 / 0.2e1 + (-Ifges(6,2) - Ifges(5,1) / 0.2e1) * t90 + (-Ifges(6,6) - Ifges(5,4)) * t91) * t192 + t265 * t400 + t355 * t419; -(-Ifges(3,2) * t328 + t219 + t254) * t329 / 0.2e1 + t424 * (-pkin(10) + t211) + (-t102 * t92 - t195 * t198 + t20 * t216 - t21 * t215 + (t182 - t103) * t93) * m(5) + (t422 + t444) * (t252 * t347 + (qJD(4) * t317 + (t260 * t395 + t317) * qJD(3)) * pkin(2)) + (t10 * t210 - t24 * t30 - t25 * t31 + t437 * t54) * m(7) + t437 * t101 + t438 * t210 + t432 * t102 + t433 * t103 + (-pkin(2) * t323 - t193) * t200 + Ifges(3,5) * t316 + (t319 - t194) * t199 + (-mrSges(3,1) * t316 + mrSges(3,2) * t315) * pkin(7) + t263 + ((t395 * t133 + t392 * t132 + (-t186 * t392 + t187 * t395) * qJD(3)) * pkin(2) - t186 * t193 - t187 * t194 - t242 * t255) * m(4) - t198 * t129 + t181 * t155 + t182 * t157 + (t431 + (-t293 / 0.2e1 + t297) * qJD(1)) * qJD(1) - t104 * t130 - t31 * t108 - t30 * t109 - t216 * t380 - t211 * t383 - t300 * t344 / 0.2e1 + t270 * t340 + (-t100 * t104 - t102 * t76 - t17 * t210 + t21 * t211 + t451 * t77) * m(6) - t271 * t342 - Ifges(3,6) * t315 - t185 * t255 + t218 * t328 / 0.2e1 + t215 * t382; -t107 * t130 - t29 * t108 - t28 * t109 - t129 * t390 - t186 * t199 + t187 * t200 - t251 * t383 + t339 * t382 - t380 * t389 + t263 + t432 * t98 + t438 * t250 + t458 * t157 + t436 * t155 + t435 * t101 + (t10 * t250 - t24 * t28 - t25 * t29 + t435 * t54) * m(7) + (-t100 * t107 - t17 * t250 + t21 * t251 + t436 * t77 - t76 * t98) * m(6) + ((-t394 * t21 + t20 * t260 + (t260 * t92 + t394 * t93) * qJD(4)) * pkin(3) - t195 * t390 - t92 * t98 - t93 * t99) * m(5) + t422 * pkin(3) * t347 + t424 * (-pkin(10) + t251); ((Ifges(6,6) / 0.2e1 + Ifges(5,4) / 0.2e1) * t282 + (Ifges(5,6) / 0.2e1 - Ifges(6,5) / 0.2e1) * t257 + t305 * mrSges(7,3) + (-Ifges(6,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 - Ifges(5,2) / 0.2e1) * t175 + t457) * t282 + (-(qJD(6) * t108 + t26) * t415 + (-t2 - t357) * mrSges(7,3)) * t261 + (-t1 * mrSges(7,3) - (-qJD(6) * t109 + t27) * t415) * t259 + (t366 / 0.2e1 + t120 / 0.2e1 + t62 / 0.2e1 - t118 / 0.2e1 - t165 / 0.2e1 - t166 / 0.2e1 + (Ifges(5,5) / 0.2e1 - Ifges(6,4) / 0.2e1) * t257 - t447) * t175 + (pkin(4) * t90 - qJ(5) * t91) * mrSges(6,1) + t432 * t93 - t433 * t92 + t349 * qJD(5) + t280 - t126 * t130 - t33 * t108 - t32 * t109 + t303 * t101 + qJ(5) * t16 - t415 * t283 + (t10 * qJ(5) - t24 * t32 - t25 * t33 + t448 * t54) * m(7) + (-pkin(4) * t21 - qJ(5) * t17 - t100 * t126 + t434 * t77 - t76 * t93) * m(6); -t383 - t349 * t257 + t304 * qJD(6) + (t130 + t304) * t282 - m(7) * (t257 * t54 + t282 * t305) + (t100 * t282 + t257 * t77 + t21) * m(6) + t450; -t54 * (mrSges(7,1) * t152 + mrSges(7,2) * t151) + (Ifges(7,1) * t151 - t374) * t409 + t63 * t408 + (Ifges(7,5) * t151 - Ifges(7,6) * t152) * t407 - t24 * t108 + t25 * t109 + (t151 * t24 + t152 * t25) * mrSges(7,3) + t343 + (-Ifges(7,2) * t152 + t146 + t64) * t410 + t426;];
tauc  = t5(:);
