% Calculate vector of inverse dynamics joint torques for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRRRR7_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:23
% EndTime: 2019-12-31 22:21:48
% DurationCPUTime: 12.27s
% Computational Cost: add. (11859->609), mult. (27447->827), div. (0->0), fcn. (19862->14), ass. (0->299)
t255 = qJ(2) + qJ(3);
t250 = qJ(4) + t255;
t235 = sin(t250);
t236 = cos(t250);
t256 = sin(qJ(5));
t392 = mrSges(6,2) * t256;
t466 = t235 * t392 + t236 * (m(6) * pkin(9) + mrSges(6,3));
t253 = qJD(2) + qJD(3);
t247 = qJD(4) + t253;
t261 = cos(qJ(5));
t258 = sin(qJ(3));
t259 = sin(qJ(2));
t263 = cos(qJ(3));
t264 = cos(qJ(2));
t202 = -t258 * t259 + t263 * t264;
t187 = t202 * qJD(1);
t203 = t258 * t264 + t259 * t263;
t188 = t203 * qJD(1);
t257 = sin(qJ(4));
t262 = cos(qJ(4));
t293 = t187 * t257 + t262 * t188;
t122 = t247 * t261 - t256 * t293;
t123 = t247 * t256 + t261 * t293;
t445 = mrSges(5,1) * t247 + mrSges(6,1) * t122 - mrSges(6,2) * t123 - mrSges(5,3) * t293;
t266 = -pkin(7) - pkin(6);
t225 = t266 * t264;
t209 = qJD(1) * t225;
t189 = t258 * t209;
t224 = t266 * t259;
t208 = qJD(1) * t224;
t196 = qJD(2) * pkin(2) + t208;
t148 = t263 * t196 + t189;
t181 = t188 * pkin(8);
t120 = t148 - t181;
t112 = pkin(3) * t253 + t120;
t192 = t263 * t209;
t149 = t196 * t258 - t192;
t402 = pkin(8) * t187;
t121 = t149 + t402;
t366 = t121 * t257;
t71 = t112 * t262 - t366;
t67 = -pkin(4) * t247 - t71;
t465 = -m(5) * t71 + m(6) * t67 - t445;
t313 = t262 * t187 - t188 * t257;
t135 = Ifges(5,4) * t313;
t251 = t264 * pkin(2);
t240 = t251 + pkin(1);
t223 = t240 * qJD(1);
t162 = -pkin(3) * t187 - t223;
t304 = mrSges(6,1) * t256 + mrSges(6,2) * t261;
t286 = t67 * t304;
t117 = Ifges(6,4) * t122;
t138 = qJD(5) - t313;
t58 = t123 * Ifges(6,1) + t138 * Ifges(6,5) + t117;
t369 = t261 * t58;
t377 = t247 * Ifges(5,5);
t397 = t71 * mrSges(5,3);
t413 = t256 / 0.2e1;
t384 = t123 * Ifges(6,4);
t57 = t122 * Ifges(6,2) + t138 * Ifges(6,6) + t384;
t93 = Ifges(5,1) * t293 + t135 + t377;
t464 = -t369 / 0.2e1 + t57 * t413 + t397 - t93 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t293 - t135 / 0.2e1 - t286 - t162 * mrSges(5,2) - t377 / 0.2e1;
t376 = t247 * Ifges(5,6);
t388 = Ifges(5,4) * t293;
t352 = t262 * t121;
t72 = t112 * t257 + t352;
t396 = t72 * mrSges(5,3);
t68 = pkin(9) * t247 + t72;
t75 = -pkin(4) * t313 - pkin(9) * t293 + t162;
t24 = t256 * t75 + t261 * t68;
t447 = t24 * mrSges(6,2);
t23 = -t256 * t68 + t261 * t75;
t448 = t23 * mrSges(6,1);
t382 = t138 * Ifges(6,3);
t383 = t123 * Ifges(6,5);
t385 = t122 * Ifges(6,6);
t56 = t382 + t383 + t385;
t92 = Ifges(5,2) * t313 + t376 + t388;
t463 = t447 + t396 + t92 / 0.2e1 - t56 / 0.2e1 + t388 / 0.2e1 - t162 * mrSges(5,1) + t376 / 0.2e1 - t448;
t442 = t236 * pkin(4) + t235 * pkin(9);
t460 = m(6) * t442;
t459 = -t236 * mrSges(5,1) + (mrSges(5,2) - mrSges(6,3)) * t235;
t248 = sin(t255);
t249 = cos(t255);
t393 = mrSges(5,2) * t236;
t458 = mrSges(4,1) * t248 + mrSges(5,1) * t235 + mrSges(4,2) * t249 + t393;
t340 = qJD(1) * qJD(2);
t212 = qJDD(1) * t264 - t259 * t340;
t371 = t261 * mrSges(6,1);
t457 = t371 - t392;
t260 = sin(qJ(1));
t265 = cos(qJ(1));
t456 = g(1) * t265 + g(2) * t260;
t450 = m(5) + m(6);
t449 = t212 / 0.2e1;
t412 = t264 / 0.2e1;
t446 = t264 * Ifges(3,2);
t164 = t258 * t224 - t263 * t225;
t315 = t249 * mrSges(4,1) - mrSges(4,2) * t248;
t441 = -t236 * t457 + t459;
t154 = t202 * t257 + t203 * t262;
t341 = qJD(5) * t261;
t280 = t202 * qJD(3);
t157 = qJD(2) * t202 + t280;
t281 = t203 * qJD(3);
t158 = -qJD(2) * t203 - t281;
t292 = t262 * t202 - t203 * t257;
t80 = qJD(4) * t292 + t157 * t262 + t158 * t257;
t288 = t154 * t341 + t256 * t80;
t342 = qJD(5) * t256;
t440 = -t23 * t341 - t24 * t342;
t348 = qJD(1) * t264;
t349 = qJD(1) * t259;
t403 = pkin(6) * t264;
t404 = pkin(6) * t259;
t439 = (qJD(2) * mrSges(3,1) - mrSges(3,3) * t349) * t403 + (-qJD(2) * mrSges(3,2) + mrSges(3,3) * t348) * t404;
t198 = t212 * pkin(6);
t213 = qJDD(1) * t259 + t264 * t340;
t199 = t213 * pkin(6);
t438 = t198 * t264 + t199 * t259;
t129 = -mrSges(5,2) * t247 + mrSges(5,3) * t313;
t86 = -mrSges(6,2) * t138 + mrSges(6,3) * t122;
t87 = mrSges(6,1) * t138 - mrSges(6,3) * t123;
t436 = -t256 * t87 + t261 * t86 + t129;
t434 = -t315 + t441;
t126 = qJD(1) * t280 + t212 * t258 + t213 * t263;
t252 = qJDD(2) + qJDD(3);
t160 = qJDD(2) * pkin(2) - pkin(7) * t213 - t199;
t161 = pkin(7) * t212 + t198;
t84 = -qJD(3) * t149 + t263 * t160 - t161 * t258;
t47 = pkin(3) * t252 - pkin(8) * t126 + t84;
t127 = -qJD(1) * t281 + t212 * t263 - t213 * t258;
t345 = qJD(3) * t263;
t346 = qJD(3) * t258;
t83 = t258 * t160 + t263 * t161 + t196 * t345 + t209 * t346;
t54 = pkin(8) * t127 + t83;
t15 = -qJD(4) * t72 - t257 * t54 + t262 * t47;
t433 = -m(3) * pkin(6) + m(4) * t266 + mrSges(2,2) - mrSges(3,3) - mrSges(4,3) - mrSges(5,3);
t63 = -qJD(4) * t293 - t126 * t257 + t127 * t262;
t246 = qJDD(4) + t252;
t62 = qJD(4) * t313 + t126 * t262 + t127 * t257;
t33 = qJD(5) * t122 + t246 * t256 + t261 * t62;
t61 = qJDD(5) - t63;
t17 = mrSges(6,1) * t61 - mrSges(6,3) * t33;
t34 = -qJD(5) * t123 + t246 * t261 - t256 * t62;
t18 = -mrSges(6,2) * t61 + mrSges(6,3) * t34;
t432 = -t256 * t17 + t261 * t18 - t87 * t341 - t86 * t342;
t343 = qJD(4) * t262;
t344 = qJD(4) * t257;
t14 = t112 * t343 - t121 * t344 + t257 * t47 + t262 * t54;
t11 = pkin(9) * t246 + t14;
t367 = qJDD(1) * pkin(1);
t179 = -pkin(2) * t212 - t367;
t102 = -pkin(3) * t127 + t179;
t19 = -pkin(4) * t63 - pkin(9) * t62 + t102;
t2 = qJD(5) * t23 + t11 * t261 + t19 * t256;
t3 = -qJD(5) * t24 - t11 * t256 + t19 * t261;
t431 = t3 * mrSges(6,1) - t2 * mrSges(6,2);
t101 = pkin(4) * t293 - pkin(9) * t313;
t222 = -mrSges(3,1) * t264 + mrSges(3,2) * t259;
t430 = m(3) * pkin(1) + m(4) * t240 + mrSges(2,1) - t222 + t315 - t459;
t298 = t23 * t261 + t24 * t256;
t398 = t3 * t256;
t273 = -qJD(5) * t298 - t398;
t399 = t2 * t261;
t429 = m(6) * (t273 + t399) + t432;
t427 = m(6) * pkin(4);
t426 = t33 / 0.2e1;
t425 = t34 / 0.2e1;
t422 = t61 / 0.2e1;
t419 = -t122 / 0.2e1;
t418 = -t123 / 0.2e1;
t417 = t123 / 0.2e1;
t416 = -t138 / 0.2e1;
t414 = t188 / 0.2e1;
t411 = pkin(2) * t259;
t410 = pkin(2) * t263;
t409 = pkin(3) * t188;
t408 = pkin(3) * t248;
t234 = pkin(3) * t249;
t407 = pkin(3) * t257;
t406 = pkin(3) * t262;
t405 = pkin(4) * t235;
t391 = mrSges(6,3) * t261;
t390 = Ifges(3,4) * t259;
t389 = Ifges(3,4) * t264;
t387 = Ifges(6,4) * t256;
t386 = Ifges(6,4) * t261;
t381 = t148 * mrSges(4,3);
t380 = t188 * mrSges(4,3);
t379 = t188 * Ifges(4,4);
t375 = t256 * mrSges(6,3);
t365 = t154 * t256;
t364 = t154 * t261;
t358 = t256 * t260;
t357 = t256 * t265;
t356 = t257 * t258;
t355 = t258 * t262;
t354 = t260 * t261;
t353 = t261 * t265;
t156 = t263 * t208 + t189;
t239 = pkin(3) + t410;
t183 = pkin(2) * t355 + t257 * t239;
t350 = t234 + t251;
t347 = qJD(2) * t259;
t338 = Ifges(6,5) * t33 + Ifges(6,6) * t34 + Ifges(6,3) * t61;
t244 = pkin(2) * t347;
t332 = t235 * t371;
t243 = pkin(2) * t349;
t324 = t369 / 0.2e1;
t323 = t234 + t442;
t322 = qJD(2) * t266;
t317 = -t342 / 0.2e1;
t143 = -pkin(3) * t158 + t244;
t316 = t340 / 0.2e1;
t155 = -t208 * t258 + t192;
t163 = t263 * t224 + t225 * t258;
t168 = -pkin(3) * t202 - t240;
t311 = t466 * t260;
t310 = t466 * t265;
t308 = mrSges(3,1) * t259 + mrSges(3,2) * t264;
t303 = Ifges(6,1) * t261 - t387;
t302 = t390 + t446;
t301 = -Ifges(6,2) * t256 + t386;
t300 = Ifges(3,5) * t264 - Ifges(3,6) * t259;
t299 = Ifges(6,5) * t261 - Ifges(6,6) * t256;
t297 = -t23 * t256 + t24 * t261;
t90 = -pkin(4) * t292 - pkin(9) * t154 + t168;
t136 = -pkin(8) * t203 + t163;
t137 = pkin(8) * t202 + t164;
t95 = t136 * t257 + t137 * t262;
t40 = -t256 * t95 + t261 * t90;
t41 = t256 * t90 + t261 * t95;
t294 = t262 * t136 - t137 * t257;
t182 = -pkin(2) * t356 + t239 * t262;
t290 = t155 - t402;
t289 = pkin(1) * t308;
t287 = t154 * t342 - t261 * t80;
t285 = t122 * t301;
t284 = t123 * t303;
t283 = t138 * t299;
t282 = t259 * (Ifges(3,1) * t264 - t390);
t210 = t259 * t322;
t211 = t264 * t322;
t106 = t263 * t210 + t258 * t211 + t224 * t345 + t225 * t346;
t85 = t101 + t409;
t214 = -t408 - t411;
t274 = m(6) * (t214 - t405) - t332;
t272 = m(6) * (-t405 - t408) - t332;
t271 = t393 + (mrSges(5,1) + t371 + t427) * t235;
t107 = -qJD(3) * t164 - t210 * t258 + t263 * t211;
t269 = -pkin(8) * t157 + t107;
t12 = -pkin(4) * t246 - t15;
t8 = t33 * Ifges(6,4) + t34 * Ifges(6,2) + t61 * Ifges(6,6);
t9 = t33 * Ifges(6,1) + t34 * Ifges(6,4) + t61 * Ifges(6,5);
t268 = -t14 * mrSges(5,2) + t2 * t391 + t261 * t8 / 0.2e1 + t9 * t413 - t12 * t457 + t15 * mrSges(5,1) + Ifges(5,3) * t246 + (Ifges(6,1) * t256 + t386) * t426 + (Ifges(6,2) * t261 + t387) * t425 + (Ifges(6,5) * t256 + Ifges(6,6) * t261) * t422 + t57 * t317 + Ifges(5,6) * t63 + Ifges(5,5) * t62 + (t324 + t286) * qJD(5) + (t285 + t284 + t283) * qJD(5) / 0.2e1;
t133 = t187 * Ifges(4,2) + t253 * Ifges(4,6) + t379;
t180 = Ifges(4,4) * t187;
t134 = t188 * Ifges(4,1) + t253 * Ifges(4,5) + t180;
t267 = t268 + t440 * mrSges(6,3) - t253 * (Ifges(4,5) * t187 - Ifges(4,6) * t188) / 0.2e1 + Ifges(4,3) * t252 + t223 * (mrSges(4,1) * t188 + mrSges(4,2) * t187) + Ifges(4,6) * t127 + Ifges(4,5) * t126 + t149 * t380 + t187 * t381 + t133 * t414 - t3 * t375 - t188 * (Ifges(4,1) * t187 - t379) / 0.2e1 - t83 * mrSges(4,2) + t84 * mrSges(4,1) - (-Ifges(4,2) * t188 + t134 + t180) * t187 / 0.2e1 + (Ifges(6,5) * t418 + Ifges(6,6) * t419 + Ifges(6,3) * t416 + t463) * t293 + (t23 * t391 + t24 * t375 + t299 * t416 + t301 * t419 + t303 * t418 + t464) * t313;
t254 = -pkin(8) + t266;
t242 = Ifges(3,4) * t348;
t238 = -pkin(4) - t406;
t207 = pkin(1) + t350;
t186 = Ifges(3,1) * t349 + Ifges(3,5) * qJD(2) + t242;
t185 = Ifges(3,6) * qJD(2) + qJD(1) * t302;
t177 = -pkin(4) - t182;
t175 = t236 * t353 + t358;
t174 = -t236 * t357 + t354;
t173 = -t236 * t354 + t357;
t172 = t236 * t358 + t353;
t167 = mrSges(4,1) * t253 - t380;
t166 = -mrSges(4,2) * t253 + mrSges(4,3) * t187;
t165 = t243 + t409;
t147 = -mrSges(4,1) * t187 + mrSges(4,2) * t188;
t128 = -t181 + t156;
t111 = -mrSges(4,2) * t252 + mrSges(4,3) * t127;
t110 = mrSges(4,1) * t252 - mrSges(4,3) * t126;
t100 = -mrSges(5,1) * t313 + mrSges(5,2) * t293;
t89 = pkin(8) * t158 + t106;
t81 = qJD(4) * t154 + t157 * t257 - t262 * t158;
t79 = t243 + t85;
t78 = t262 * t128 + t257 * t290;
t74 = t120 * t262 - t366;
t73 = t120 * t257 + t352;
t50 = -mrSges(5,2) * t246 + mrSges(5,3) * t63;
t49 = mrSges(5,1) * t246 - mrSges(5,3) * t62;
t36 = t101 * t256 + t261 * t71;
t35 = t101 * t261 - t256 * t71;
t30 = t256 * t79 + t261 * t78;
t29 = -t256 * t78 + t261 * t79;
t28 = t256 * t85 + t261 * t74;
t27 = -t256 * t74 + t261 * t85;
t22 = pkin(4) * t81 - pkin(9) * t80 + t143;
t20 = qJD(4) * t294 + t257 * t269 + t262 * t89;
t16 = -mrSges(6,1) * t34 + mrSges(6,2) * t33;
t5 = -qJD(5) * t41 - t20 * t256 + t22 * t261;
t4 = qJD(5) * t40 + t20 * t261 + t22 * t256;
t1 = [-t81 * t447 + (t212 * t403 + t213 * t404 + t438) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t438) + (t186 * t412 + t300 * qJD(2) / 0.2e1 - t439) * qJD(2) - t288 * t57 / 0.2e1 + (-mrSges(4,1) * t179 + mrSges(4,3) * t83 + Ifges(4,4) * t126 + Ifges(4,2) * t127 + Ifges(4,6) * t252) * t202 + t213 * t389 / 0.2e1 + t253 * (Ifges(4,5) * t157 + Ifges(4,6) * t158) / 0.2e1 + t247 * (Ifges(5,5) * t80 - Ifges(5,6) * t81) / 0.2e1 - t240 * (-mrSges(4,1) * t127 + mrSges(4,2) * t126) - (-Ifges(5,4) * t62 - Ifges(5,2) * t63 - Ifges(5,6) * t246 + t102 * mrSges(5,1) - t14 * mrSges(5,3) + t338 / 0.2e1 + Ifges(6,3) * t422 + Ifges(6,6) * t425 + Ifges(6,5) * t426 + t431) * t292 + (-t173 * mrSges(6,1) - t172 * mrSges(6,2) + (t254 * t450 + t433) * t265 + (m(5) * t207 - m(6) * (-t207 - t442) + t430) * t260) * g(1) + m(4) * (t106 * t149 + t107 * t148 + t163 * t84 + t164 * t83 - t179 * t240 - t223 * t244) + (t102 * mrSges(5,2) - t15 * mrSges(5,3) + Ifges(5,1) * t62 + Ifges(5,4) * t63 + Ifges(5,5) * t246 + t12 * t304 + t299 * t422 + t301 * t425 + t303 * t426 + t317 * t58) * t154 - (-m(5) * t15 + m(6) * t12 + t16 - t49) * t294 + (Ifges(3,4) * t213 + Ifges(3,2) * t212) * t412 + m(5) * (t102 * t168 + t14 * t95 + t143 * t162 + t20 * t72) + m(6) * (t2 * t41 + t23 * t5 + t24 * t4 + t3 * t40) + (t264 * t389 + t282) * t316 + t81 * t448 + t302 * t449 + t313 * (Ifges(5,4) * t80 - Ifges(5,2) * t81) / 0.2e1 + t293 * (Ifges(5,1) * t80 - Ifges(5,4) * t81) / 0.2e1 + t138 * (-Ifges(6,5) * t287 - Ifges(6,6) * t288 + Ifges(6,3) * t81) / 0.2e1 + t122 * (-Ifges(6,4) * t287 - Ifges(6,2) * t288 + Ifges(6,6) * t81) / 0.2e1 + t67 * (mrSges(6,1) * t288 - mrSges(6,2) * t287) + t147 * t244 + (mrSges(4,2) * t179 - mrSges(4,3) * t84 + Ifges(4,1) * t126 + Ifges(4,4) * t127 + Ifges(4,5) * t252) * t203 - t223 * (-mrSges(4,1) * t158 + mrSges(4,2) * t157) - pkin(1) * (-mrSges(3,1) * t212 + mrSges(3,2) * t213) + t187 * (Ifges(4,4) * t157 + Ifges(4,2) * t158) / 0.2e1 + (-t175 * mrSges(6,1) - t174 * mrSges(6,2) - t450 * (t265 * t207 - t254 * t260) + t433 * t260 + (-t430 - t460) * t265) * g(2) + t168 * (-mrSges(5,1) * t63 + mrSges(5,2) * t62) + t163 * t110 + t164 * t111 + t106 * t166 + t107 * t167 + (-t2 * t365 + t23 * t287 - t24 * t288 - t3 * t364) * mrSges(6,3) + t157 * t134 / 0.2e1 + t158 * t133 / 0.2e1 + t162 * (mrSges(5,1) * t81 + mrSges(5,2) * t80) + Ifges(2,3) * qJDD(1) + t143 * t100 + t20 * t129 + t149 * t158 * mrSges(4,3) + (-mrSges(3,1) * t404 - mrSges(3,2) * t403 + 0.2e1 * Ifges(3,6) * t412) * qJDD(2) + (Ifges(3,1) * t213 + Ifges(3,4) * t449 + Ifges(3,5) * qJDD(2) - t316 * t446) * t259 - t289 * t340 + t80 * t324 - t185 * t347 / 0.2e1 + (Ifges(4,1) * t157 + Ifges(4,4) * t158) * t414 + (-Ifges(6,1) * t287 - Ifges(6,4) * t288 + Ifges(6,5) * t81) * t417 + t9 * t364 / 0.2e1 - t8 * t365 / 0.2e1 - t222 * t367 - t157 * t381 - t81 * t396 - t80 * t397 + t465 * (qJD(4) * t95 + t257 * t89 - t262 * t269) + t40 * t17 + t41 * t18 + t81 * t56 / 0.2e1 + t4 * t86 + t5 * t87 - t81 * t92 / 0.2e1 + t80 * t93 / 0.2e1 + t95 * t50; (t222 - m(6) * (t251 + t323) - m(5) * t350 - m(4) * t251 + t434) * g(3) + t267 + (t439 + (-t282 / 0.2e1 + t289) * qJD(1)) * qJD(1) - m(4) * (t148 * t155 + t149 * t156 - t223 * t243) + (m(5) * t72 + m(6) * t297 + t436) * (t239 * t343 + (-t258 * t344 + (t262 * t263 - t356) * qJD(3)) * pkin(2)) - (-Ifges(3,2) * t349 + t186 + t242) * t348 / 0.2e1 + t429 * (pkin(9) + t183) - g(1) * (t265 * t274 + t310) + t456 * (m(4) * t411 - m(5) * t214 + t308 + t458) + (t166 * t345 + m(4) * (t258 * t83 + t263 * t84 + (-t148 * t258 + t149 * t263) * qJD(3)) - t167 * t346 + t258 * t111) * pkin(2) + Ifges(3,6) * t212 + Ifges(3,5) * t213 - t198 * mrSges(3,2) - t199 * mrSges(3,1) + t177 * t16 + t182 * t49 + t183 * t50 - t165 * t100 - t156 * t166 - t155 * t167 - t78 * t129 + (t12 * t177 - t23 * t29 - t24 * t30) * m(6) + (t14 * t183 + t15 * t182 - t162 * t165 - t72 * t78) * m(5) - t147 * t243 - t300 * t340 / 0.2e1 + t185 * t349 / 0.2e1 + t110 * t410 - g(2) * (t260 * t274 + t311) + Ifges(3,3) * qJDD(2) + t465 * (-t128 * t257 + t262 * t290 + t239 * t344 + (t258 * t343 + (t257 * t263 + t355) * qJD(3)) * pkin(2)) - t30 * t86 - t29 * t87; t267 + t238 * t16 - g(1) * (t265 * t272 + t310) - t148 * t166 + t149 * t167 - t74 * t129 + t49 * t406 + t50 * t407 - g(2) * (t260 * t272 + t311) - t100 * t409 - t28 * t86 - t27 * t87 + (t12 * t238 + (t257 * t67 + t262 * t297) * qJD(4) * pkin(3) - t23 * t27 - t24 * t28 - t67 * t73) * m(6) + ((t14 * t257 + t15 * t262 + (-t257 * t71 + t262 * t72) * qJD(4)) * pkin(3) - t162 * t409 + t71 * t73 - t72 * t74) * m(5) + t436 * pkin(3) * t343 + (-m(5) * t234 - m(6) * t323 + t434) * g(3) + t429 * (pkin(9) + t407) + t445 * (-pkin(3) * t344 + t73) + (m(5) * t408 + t458) * t456; t445 * t72 + t268 - m(6) * (t23 * t35 + t24 * t36 + t67 * t72) + (-t383 / 0.2e1 - t382 / 0.2e1 - t385 / 0.2e1 + t463) * t293 + (t441 - t460) * g(3) - t12 * t427 - t71 * t129 + (m(6) * (-t398 + t399 + t440) + t432) * pkin(9) + (t260 * t271 - t311) * g(2) + (t265 * t271 - t310) * g(1) + t273 * mrSges(6,3) + (-t284 / 0.2e1 - t283 / 0.2e1 - t285 / 0.2e1 + t298 * mrSges(6,3) + t464) * t313 - pkin(4) * t16 - t36 * t86 - t35 * t87; -t67 * (mrSges(6,1) * t123 + mrSges(6,2) * t122) + (Ifges(6,1) * t122 - t384) * t418 + t57 * t417 + (Ifges(6,5) * t122 - Ifges(6,6) * t123) * t416 - t23 * t86 + t24 * t87 - g(1) * (mrSges(6,1) * t174 - mrSges(6,2) * t175) - g(2) * (-mrSges(6,1) * t172 + mrSges(6,2) * t173) + g(3) * t304 * t235 + (t122 * t23 + t123 * t24) * mrSges(6,3) + t338 + (-Ifges(6,2) * t123 + t117 + t58) * t419 + t431;];
tau = t1;
