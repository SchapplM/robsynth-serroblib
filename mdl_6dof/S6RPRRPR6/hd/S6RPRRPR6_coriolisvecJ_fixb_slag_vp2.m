% Calculate vector of centrifugal and Coriolis load on the joints for
% S6RPRRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 05:18
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S6RPRRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:15:36
% EndTime: 2019-03-09 05:16:01
% DurationCPUTime: 13.58s
% Computational Cost: add. (16114->643), mult. (42085->883), div. (0->0), fcn. (32757->10), ass. (0->282)
t265 = sin(pkin(10));
t267 = cos(pkin(10));
t270 = sin(qJ(3));
t273 = cos(qJ(3));
t404 = -t265 * t270 + t273 * t267;
t238 = t404 * qJD(1);
t264 = sin(pkin(11));
t266 = cos(pkin(11));
t269 = sin(qJ(4));
t272 = cos(qJ(4));
t246 = t264 * t272 + t266 * t269;
t168 = t246 * t238;
t237 = t246 * qJD(4);
t312 = t168 - t237;
t283 = t264 * t269 - t266 * t272;
t169 = t283 * t238;
t239 = t283 * qJD(4);
t311 = t169 - t239;
t268 = sin(qJ(6));
t271 = cos(qJ(6));
t228 = qJD(4) - t238;
t247 = t265 * t273 + t267 * t270;
t240 = t247 * qJD(1);
t206 = qJD(3) * t272 - t240 * t269;
t207 = qJD(3) * t269 + t240 * t272;
t150 = t206 * t264 + t207 * t266;
t407 = pkin(9) * t150;
t302 = -pkin(2) * t267 - pkin(1);
t253 = qJD(1) * t302 + qJD(2);
t165 = -pkin(3) * t238 - pkin(8) * t240 + t253;
t342 = pkin(7) + qJ(2);
t254 = t342 * t265;
t251 = qJD(1) * t254;
t255 = t342 * t267;
t252 = qJD(1) * t255;
t198 = -t270 * t251 + t273 * t252;
t191 = qJD(3) * pkin(8) + t198;
t121 = t272 * t165 - t191 * t269;
t100 = -qJ(5) * t207 + t121;
t82 = pkin(4) * t228 + t100;
t122 = t165 * t269 + t191 * t272;
t101 = qJ(5) * t206 + t122;
t93 = t264 * t101;
t48 = t266 * t82 - t93;
t32 = pkin(5) * t228 - t407 + t48;
t297 = t266 * t206 - t207 * t264;
t401 = pkin(9) * t297;
t317 = t266 * t101;
t49 = t264 * t82 + t317;
t33 = t49 + t401;
t11 = -t268 * t33 + t271 * t32;
t12 = t268 * t32 + t271 * t33;
t242 = t247 * qJD(3);
t225 = qJD(1) * t242;
t241 = t404 * qJD(3);
t224 = qJD(1) * t241;
t159 = qJD(4) * t206 + t224 * t272;
t160 = -qJD(4) * t207 - t224 * t269;
t106 = -t159 * t264 + t160 * t266;
t107 = t159 * t266 + t160 * t264;
t403 = -t150 * t268 + t271 * t297;
t30 = qJD(6) * t403 + t106 * t268 + t107 * t271;
t92 = t150 * t271 + t268 * t297;
t31 = -qJD(6) * t92 + t106 * t271 - t107 * t268;
t307 = Ifges(7,5) * t30 + Ifges(7,6) * t31 + Ifges(7,3) * t225;
t352 = Ifges(7,4) * t92;
t223 = qJD(6) + t228;
t360 = -t223 / 0.2e1;
t374 = -t92 / 0.2e1;
t376 = -t403 / 0.2e1;
t197 = -t251 * t273 - t270 * t252;
t278 = t404 * qJD(2);
t152 = qJD(1) * t278 + qJD(3) * t197;
t177 = pkin(3) * t225 - pkin(8) * t224;
t71 = -qJD(4) * t122 - t152 * t269 + t272 * t177;
t38 = pkin(4) * t225 - qJ(5) * t159 - qJD(5) * t207 + t71;
t308 = qJD(4) * t272;
t309 = qJD(4) * t269;
t70 = t272 * t152 + t165 * t308 + t269 * t177 - t191 * t309;
t41 = qJ(5) * t160 + qJD(5) * t206 + t70;
t16 = t264 * t38 + t266 * t41;
t10 = pkin(9) * t106 + t16;
t15 = -t264 * t41 + t266 * t38;
t8 = pkin(5) * t225 - pkin(9) * t107 + t15;
t2 = qJD(6) * t11 + t10 * t271 + t268 * t8;
t3 = -qJD(6) * t12 - t10 * t268 + t271 * t8;
t400 = t3 * mrSges(7,1) - t2 * mrSges(7,2);
t87 = Ifges(7,4) * t403;
t46 = Ifges(7,1) * t92 + Ifges(7,5) * t223 + t87;
t190 = -qJD(3) * pkin(3) - t197;
t144 = -pkin(4) * t206 + qJD(5) + t190;
t88 = -pkin(5) * t297 + t144;
t421 = t307 + t400 + (Ifges(7,5) * t403 - Ifges(7,6) * t92) * t360 + (t11 * t403 + t12 * t92) * mrSges(7,3) + (-Ifges(7,2) * t92 + t46 + t87) * t376 - t88 * (mrSges(7,1) * t92 + mrSges(7,2) * t403) + (Ifges(7,1) * t403 - t352) * t374;
t341 = -qJ(5) - pkin(8);
t300 = qJD(4) * t341;
t232 = qJD(5) * t272 + t269 * t300;
t233 = -qJD(5) * t269 + t272 * t300;
t173 = -t232 * t264 + t266 * t233;
t192 = pkin(3) * t240 - pkin(8) * t238;
t133 = t272 * t192 - t197 * t269;
t325 = qJ(5) * t272;
t102 = pkin(4) * t240 - t238 * t325 + t133;
t134 = t269 * t192 + t272 * t197;
t322 = t238 * t269;
t117 = -qJ(5) * t322 + t134;
t61 = t266 * t102 - t117 * t264;
t420 = t173 - t61;
t174 = t266 * t232 + t264 * t233;
t62 = t264 * t102 + t266 * t117;
t419 = t174 - t62;
t45 = Ifges(7,2) * t403 + Ifges(7,6) * t223 + t352;
t417 = t45 / 0.2e1;
t416 = Ifges(5,3) + Ifges(6,3);
t411 = -pkin(5) * t240 - pkin(9) * t311 + t420;
t410 = pkin(9) * t312 + t419;
t335 = t207 * Ifges(5,4);
t138 = t206 * Ifges(5,2) + t228 * Ifges(5,6) + t335;
t205 = Ifges(5,4) * t206;
t139 = t207 * Ifges(5,1) + t228 * Ifges(5,5) + t205;
t286 = t121 * t272 + t122 * t269;
t289 = Ifges(5,5) * t272 - Ifges(5,6) * t269;
t339 = Ifges(5,4) * t272;
t291 = -Ifges(5,2) * t269 + t339;
t340 = Ifges(5,4) * t269;
t293 = Ifges(5,1) * t272 - t340;
t294 = mrSges(5,1) * t269 + mrSges(5,2) * t272;
t353 = t272 / 0.2e1;
t354 = -t269 / 0.2e1;
t356 = t228 / 0.2e1;
t361 = t207 / 0.2e1;
t409 = t289 * t356 + t206 * t291 / 0.2e1 + t293 * t361 + t190 * t294 + t138 * t354 + t139 * t353 - t286 * mrSges(5,3);
t408 = -t238 * Ifges(4,2) / 0.2e1;
t406 = t150 * Ifges(6,4);
t227 = Ifges(4,4) * t238;
t387 = t227 / 0.2e1 + t240 * Ifges(4,1) / 0.2e1;
t405 = t253 * mrSges(4,2) + Ifges(4,5) * qJD(3) + t387 + t409;
t380 = t30 / 0.2e1;
t379 = t31 / 0.2e1;
t80 = Ifges(6,2) * t297 + t228 * Ifges(6,6) + t406;
t402 = t80 / 0.2e1;
t371 = t106 / 0.2e1;
t370 = t107 / 0.2e1;
t358 = t225 / 0.2e1;
t357 = -t228 / 0.2e1;
t369 = -t297 / 0.2e1;
t256 = t341 * t269;
t257 = t341 * t272;
t203 = t266 * t256 + t257 * t264;
t175 = -pkin(9) * t246 + t203;
t204 = t264 * t256 - t266 * t257;
t176 = -pkin(9) * t283 + t204;
t119 = t175 * t271 - t176 * t268;
t399 = qJD(6) * t119 + t268 * t411 + t271 * t410;
t120 = t175 * t268 + t176 * t271;
t398 = -qJD(6) * t120 - t268 * t410 + t271 * t411;
t397 = Ifges(6,4) * t297;
t260 = pkin(4) * t266 + pkin(5);
t350 = pkin(4) * t264;
t230 = t260 * t271 - t268 * t350;
t53 = -t100 * t264 - t317;
t34 = t53 - t401;
t54 = t266 * t100 - t93;
t35 = t54 - t407;
t396 = t230 * qJD(6) - t268 * t34 - t271 * t35;
t231 = t260 * t268 + t271 * t350;
t395 = -t231 * qJD(6) + t268 * t35 - t271 * t34;
t116 = -t168 * t268 - t169 * t271;
t195 = -t246 * t268 - t271 * t283;
t135 = qJD(6) * t195 - t237 * t268 - t239 * t271;
t315 = t135 - t116;
t115 = -t168 * t271 + t169 * t268;
t196 = t246 * t271 - t268 * t283;
t136 = -qJD(6) * t196 - t237 * t271 + t239 * t268;
t314 = t136 - t115;
t392 = Ifges(5,5) * t159 + Ifges(5,6) * t160;
t194 = -pkin(3) * t404 - pkin(8) * t247 + t302;
t202 = -t254 * t270 + t255 * t273;
t199 = t272 * t202;
t143 = t269 * t194 + t199;
t155 = pkin(4) * t322 + t198;
t391 = pkin(4) * t309 - pkin(5) * t312 - t155;
t390 = -t273 * t254 - t255 * t270;
t389 = Ifges(6,5) * t107 + Ifges(6,6) * t106 + t225 * t416 + t392;
t388 = -t269 * t71 + t272 * t70;
t386 = (m(3) * qJ(2) + mrSges(3,3)) * (t265 ^ 2 + t267 ^ 2);
t381 = Ifges(5,3) / 0.2e1;
t385 = (t381 + Ifges(6,3) / 0.2e1) * t228 - t12 * mrSges(7,2) - t122 * mrSges(5,2) - t49 * mrSges(6,2) + t207 * Ifges(5,5) + t206 * Ifges(5,6) - Ifges(4,6) * qJD(3) - t240 * Ifges(4,4) + t408 + t223 * Ifges(7,3) + t92 * Ifges(7,5) + t403 * Ifges(7,6) + t150 * Ifges(6,5) + t297 * Ifges(6,6) + t11 * mrSges(7,1) + t121 * mrSges(5,1) + t253 * mrSges(4,1) + t48 * mrSges(6,1) - t416 * t357;
t384 = t71 * mrSges(5,1) + t15 * mrSges(6,1) - t70 * mrSges(5,2) - t16 * mrSges(6,2);
t383 = Ifges(7,4) * t380 + Ifges(7,2) * t379 + Ifges(7,6) * t358;
t382 = Ifges(7,1) * t380 + Ifges(7,4) * t379 + Ifges(7,5) * t358;
t378 = Ifges(6,4) * t370 + Ifges(6,2) * t371 + Ifges(6,6) * t358;
t377 = Ifges(6,1) * t370 + Ifges(6,4) * t371 + Ifges(6,5) * t358;
t375 = t403 / 0.2e1;
t373 = t92 / 0.2e1;
t368 = t297 / 0.2e1;
t367 = -t150 / 0.2e1;
t366 = t150 / 0.2e1;
t365 = t159 / 0.2e1;
t364 = t160 / 0.2e1;
t363 = -t206 / 0.2e1;
t362 = -t207 / 0.2e1;
t359 = t223 / 0.2e1;
t351 = pkin(4) * t207;
t349 = pkin(4) * t269;
t282 = -qJ(5) * t241 - qJD(5) * t247;
t166 = qJD(3) * t390 + t278;
t193 = pkin(3) * t242 - pkin(8) * t241;
t298 = -t166 * t269 + t272 * t193;
t59 = pkin(4) * t242 + t282 * t272 + (-t199 + (qJ(5) * t247 - t194) * t269) * qJD(4) + t298;
t301 = t247 * t308;
t303 = t272 * t166 + t269 * t193 + t194 * t308;
t63 = -qJ(5) * t301 + (-qJD(4) * t202 + t282) * t269 + t303;
t24 = t264 * t59 + t266 * t63;
t279 = t247 * qJD(2);
t153 = qJD(1) * t279 + qJD(3) * t198;
t324 = t153 * t390;
t321 = t247 * t269;
t142 = t272 * t194 - t202 * t269;
t112 = -pkin(4) * t404 - t247 * t325 + t142;
t123 = -qJ(5) * t321 + t143;
t69 = t264 * t112 + t266 * t123;
t313 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t206 - mrSges(5,2) * t207 - mrSges(4,3) * t240;
t261 = -pkin(4) * t272 - pkin(3);
t9 = -t31 * mrSges(7,1) + t30 * mrSges(7,2);
t23 = -t264 * t63 + t266 * t59;
t299 = t225 * mrSges(4,1) + t224 * mrSges(4,2);
t66 = -t106 * mrSges(6,1) + t107 * mrSges(6,2);
t68 = t266 * t112 - t123 * t264;
t170 = pkin(4) * t321 - t390;
t295 = mrSges(5,1) * t272 - mrSges(5,2) * t269;
t292 = Ifges(5,1) * t269 + t339;
t290 = Ifges(5,2) * t272 + t340;
t288 = Ifges(5,5) * t269 + Ifges(5,6) * t272;
t180 = t283 * t247;
t50 = -pkin(5) * t404 + pkin(9) * t180 + t68;
t179 = t246 * t247;
t52 = -pkin(9) * t179 + t69;
t21 = -t268 * t52 + t271 * t50;
t22 = t268 * t50 + t271 * t52;
t287 = -t269 * t70 - t272 * t71;
t285 = t121 * t269 - t122 * t272;
t163 = -mrSges(5,2) * t228 + mrSges(5,3) * t206;
t164 = mrSges(5,1) * t228 - mrSges(5,3) * t207;
t284 = t163 * t272 - t164 * t269;
t126 = -t179 * t271 + t180 * t268;
t127 = -t179 * t268 - t180 * t271;
t167 = qJD(3) * t202 + t279;
t130 = pkin(4) * t301 + t241 * t349 + t167;
t110 = -pkin(4) * t160 + t153;
t211 = pkin(5) * t283 + t261;
t209 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t238;
t141 = -mrSges(5,2) * t225 + mrSges(5,3) * t160;
t140 = mrSges(5,1) * t225 - mrSges(5,3) * t159;
t132 = mrSges(6,1) * t228 - mrSges(6,3) * t150;
t131 = -mrSges(6,2) * t228 + mrSges(6,3) * t297;
t129 = pkin(5) * t150 + t351;
t128 = pkin(5) * t179 + t170;
t125 = -t237 * t247 - t241 * t283;
t124 = t239 * t247 - t241 * t246;
t111 = -mrSges(5,1) * t160 + mrSges(5,2) * t159;
t97 = -mrSges(6,1) * t297 + mrSges(6,2) * t150;
t86 = t159 * Ifges(5,1) + t160 * Ifges(5,4) + t225 * Ifges(5,5);
t85 = t159 * Ifges(5,4) + t160 * Ifges(5,2) + t225 * Ifges(5,6);
t84 = mrSges(6,1) * t225 - mrSges(6,3) * t107;
t83 = -mrSges(6,2) * t225 + mrSges(6,3) * t106;
t81 = t150 * Ifges(6,1) + t228 * Ifges(6,5) + t397;
t76 = mrSges(7,1) * t223 - mrSges(7,3) * t92;
t75 = -mrSges(7,2) * t223 + mrSges(7,3) * t403;
t74 = -qJD(4) * t143 + t298;
t73 = -t202 * t309 + t303;
t72 = -pkin(5) * t124 + t130;
t67 = -pkin(5) * t106 + t110;
t51 = -mrSges(7,1) * t403 + mrSges(7,2) * t92;
t43 = -qJD(6) * t127 + t124 * t271 - t125 * t268;
t42 = qJD(6) * t126 + t124 * t268 + t125 * t271;
t26 = -mrSges(7,2) * t225 + mrSges(7,3) * t31;
t25 = mrSges(7,1) * t225 - mrSges(7,3) * t30;
t20 = pkin(9) * t124 + t24;
t17 = pkin(5) * t242 - pkin(9) * t125 + t23;
t5 = -qJD(6) * t22 + t17 * t271 - t20 * t268;
t4 = qJD(6) * t21 + t17 * t268 + t20 * t271;
t1 = [-(mrSges(4,3) * t224 + t111) * t390 - (-mrSges(4,3) * t152 - Ifges(4,4) * t224 + Ifges(6,5) * t370 + Ifges(7,5) * t380 + Ifges(6,6) * t371 + Ifges(7,6) * t379 + (Ifges(6,3) + Ifges(7,3)) * t358 + t384 + t400) * t404 + (Ifges(4,1) * t224 - Ifges(4,4) * t225 + t293 * t365 + t291 * t364 + t85 * t354 + t86 * t353 + (mrSges(4,3) + t294) * t153 + t287 * mrSges(5,3) + (t190 * t295 + t288 * t357 + t290 * t363 + t292 * t362 + t139 * t354 - t272 * t138 / 0.2e1 + t285 * mrSges(5,3)) * qJD(4)) * t247 + (Ifges(7,4) * t127 + Ifges(7,2) * t126) * t379 + m(5) * (t121 * t74 + t122 * t73 + t142 * t71 + t143 * t70 + t167 * t190 - t324) + m(4) * (t152 * t202 + t166 * t198 - t167 * t197 - t324) + (Ifges(7,4) * t42 + Ifges(7,2) * t43) * t375 - t180 * t377 - t179 * t378 + t127 * t382 + t126 * t383 + (-t197 * t241 - t198 * t242 - t202 * t225) * mrSges(4,3) + (-Ifges(6,1) * t180 - Ifges(6,4) * t179) * t370 + (t124 * t49 - t125 * t48 + t15 * t180 - t16 * t179) * mrSges(6,3) + (-Ifges(6,4) * t180 - Ifges(6,2) * t179) * t371 + (-Ifges(6,5) * t180 + Ifges(7,5) * t127 - Ifges(6,6) * t179 + Ifges(7,6) * t126 + t289 * t247) * t358 + t110 * (mrSges(6,1) * t179 - mrSges(6,2) * t180) + (Ifges(7,5) * t42 + Ifges(7,6) * t43) * t359 + (Ifges(6,1) * t125 + Ifges(6,4) * t124) * t366 + (Ifges(6,4) * t125 + Ifges(6,2) * t124) * t368 + (Ifges(7,1) * t42 + Ifges(7,4) * t43) * t373 + (Ifges(6,5) * t125 + Ifges(6,6) * t124) * t356 - (t381 + Ifges(4,2)) * t225 * t404 + (-t11 * t42 + t12 * t43 + t126 * t2 - t127 * t3) * mrSges(7,3) + (Ifges(7,1) * t127 + Ifges(7,4) * t126) * t380 + 0.2e1 * t386 * qJD(2) * qJD(1) + m(7) * (t11 * t5 + t12 * t4 + t128 * t67 + t2 * t22 + t21 * t3 + t72 * t88) + m(6) * (t110 * t170 + t130 * t144 + t15 * t68 + t16 * t69 + t23 * t48 + t24 * t49) - (t307 + t389 + t392) * t404 / 0.2e1 + (t387 + t405) * t241 + t166 * t209 + t170 * t66 + t73 * t163 + t74 * t164 + t142 * t140 + t143 * t141 + t144 * (-mrSges(6,1) * t124 + mrSges(6,2) * t125) + t128 * t9 + t130 * t97 + t24 * t131 + t23 * t132 + t67 * (-mrSges(7,1) * t126 + mrSges(7,2) * t127) + t125 * t81 / 0.2e1 + t88 * (-mrSges(7,1) * t43 + mrSges(7,2) * t42) + t68 * t84 + t69 * t83 + t72 * t51 + t4 * t75 + t5 * t76 + t42 * t46 / 0.2e1 + t21 * t25 + t22 * t26 + t124 * t402 + (t408 + t385) * t242 + t302 * t299 + t43 * t417 - t313 * t167; (-t209 - t284) * t238 + t284 * qJD(4) + t312 * t132 + t311 * t131 + t269 * t141 - m(4) * (-t197 * t240 + t198 * t238) + t314 * t76 + t315 * t75 - t283 * t84 + t195 * t25 + t196 * t26 + t246 * t83 + t299 + t272 * t140 + (-t97 - t51 + t313) * t240 - t386 * qJD(1) ^ 2 + (t11 * t314 + t12 * t315 + t195 * t3 + t196 * t2 - t240 * t88) * m(7) + (-t144 * t240 - t15 * t283 + t16 * t246 + t311 * t49 + t312 * t48) * m(6) + (-t190 * t240 - t228 * t285 - t287) * m(5); -t283 * t378 + (Ifges(6,5) * t246 + Ifges(7,5) * t196 - Ifges(6,6) * t283 + Ifges(7,6) * t195 + t288) * t358 + (Ifges(6,1) * t246 - Ifges(6,4) * t283) * t370 + (Ifges(6,4) * t246 - Ifges(6,2) * t283) * t371 + t110 * (mrSges(6,1) * t283 + mrSges(6,2) * t246) + (-t15 * t246 - t16 * t283 - t311 * t48 + t312 * t49) * mrSges(6,3) + (t135 / 0.2e1 - t116 / 0.2e1) * t46 + (t136 / 0.2e1 - t115 / 0.2e1) * t45 + (Ifges(7,1) * t116 + Ifges(7,4) * t115) * t374 + (Ifges(7,4) * t135 + Ifges(7,2) * t136) * t375 + (Ifges(7,4) * t116 + Ifges(7,2) * t115) * t376 + t246 * t377 + (Ifges(7,4) * t196 + Ifges(7,2) * t195) * t379 + (Ifges(7,1) * t196 + Ifges(7,4) * t195) * t380 + t196 * t382 + t195 * t383 + m(6) * (t110 * t261 + t15 * t203 + t16 * t204 + t173 * t48 + t174 * t49) + t419 * t131 + t420 * t132 + (-Ifges(6,1) * t169 - Ifges(6,4) * t168) * t367 + (-Ifges(6,4) * t169 - Ifges(6,2) * t168) * t369 + (-Ifges(6,5) * t169 - Ifges(6,6) * t168) * t357 + (-pkin(3) * t153 - t121 * t133 - t122 * t134 - t190 * t198) * m(5) + (-mrSges(4,1) - t295) * t153 + (Ifges(7,5) * t135 + Ifges(7,6) * t136) * t359 + (Ifges(7,5) * t116 + Ifges(7,6) * t115) * t360 + t290 * t364 + t292 * t365 + (Ifges(7,1) * t135 + Ifges(7,4) * t136) * t373 + t85 * t353 + (-t237 / 0.2e1 + t168 / 0.2e1) * t80 + (-t239 / 0.2e1 + t169 / 0.2e1) * t81 + (-Ifges(6,1) * t239 - Ifges(6,4) * t237) * t366 + (-Ifges(6,4) * t239 - Ifges(6,2) * t237) * t368 + (-Ifges(6,5) * t239 - Ifges(6,6) * t237) * t356 + t261 * t66 + t269 * t86 / 0.2e1 + ((m(6) * t144 + t97) * t349 + t409) * qJD(4) + (-t140 * t269 + t141 * t272 + m(5) * t388 + (-m(5) * t286 - t269 * t163 - t272 * t164) * qJD(4)) * pkin(8) + t388 * mrSges(5,3) + t391 * t51 + (t198 * mrSges(4,3) - t385) * t240 - m(6) * (t144 * t155 + t48 * t61 + t49 * t62) + t398 * t76 + t399 * t75 + (t11 * t398 + t119 * t3 + t12 * t399 + t120 * t2 + t211 * t67 + t391 * t88) * m(7) + (-t227 / 0.2e1 + t197 * mrSges(4,3) + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t240 - t405) * t238 + Ifges(4,5) * t224 - Ifges(4,6) * t225 - t197 * t209 + t211 * t9 + t203 * t84 + t204 * t83 + t67 * (-mrSges(7,1) * t195 + mrSges(7,2) * t196) - t134 * t163 - t133 * t164 - t155 * t97 - t152 * mrSges(4,2) + t119 * t25 + t120 * t26 - pkin(3) * t111 + (-mrSges(6,1) * t312 + mrSges(6,2) * t311) * t144 + t313 * t198 + (-mrSges(7,1) * t314 + mrSges(7,2) * t315) * t88 + (-t11 * t315 + t12 * t314 + t195 * t2 - t196 * t3) * mrSges(7,3); (Ifges(5,5) * t206 + Ifges(6,5) * t297 - Ifges(5,6) * t207 - Ifges(6,6) * t150) * t357 + (t150 * t49 + t297 * t48) * mrSges(6,3) - t144 * (mrSges(6,1) * t150 + mrSges(6,2) * t297) + (-Ifges(6,2) * t150 + t397 + t81) * t369 + t150 * t402 + (-Ifges(5,2) * t207 + t139 + t205) * t363 + t389 + t138 * t361 + (Ifges(5,1) * t206 - t335) * t362 + t395 * t76 + (t11 * t395 + t12 * t396 - t129 * t88 + t2 * t231 + t230 * t3) * m(7) + t396 * t75 + t421 + t384 + ((t15 * t266 + t16 * t264) * pkin(4) - t144 * t351 - t48 * t53 - t49 * t54) * m(6) + (Ifges(6,1) * t297 - t406) * t367 + t92 * t417 + (-t207 * t97 + t264 * t83 + t266 * t84) * pkin(4) + t230 * t25 + t231 * t26 + (t121 * t206 + t122 * t207) * mrSges(5,3) - t190 * (mrSges(5,1) * t207 + mrSges(5,2) * t206) - t121 * t163 + t122 * t164 - t129 * t51 - t54 * t131 - t53 * t132; -t297 * t131 + t150 * t132 - t403 * t75 + t92 * t76 + t66 + t9 + (t11 * t92 - t12 * t403 + t67) * m(7) + (t150 * t48 - t297 * t49 + t110) * m(6); -t11 * t75 + t12 * t76 + t45 * t373 + t421;];
tauc  = t1(:);
