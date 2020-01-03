% Calculate vector of inverse dynamics joint torques for
% S5RRPRP9
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 20:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP9_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP9_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:05:42
% EndTime: 2019-12-31 20:06:14
% DurationCPUTime: 16.68s
% Computational Cost: add. (4801->554), mult. (10866->711), div. (0->0), fcn. (7397->10), ass. (0->269)
t379 = mrSges(5,1) + mrSges(6,1);
t378 = mrSges(5,2) - mrSges(6,3);
t376 = -Ifges(5,4) + Ifges(6,5);
t401 = t376 + Ifges(6,5);
t212 = cos(qJ(2));
t281 = qJD(1) * t212;
t188 = qJD(4) - t281;
t322 = -t188 / 0.2e1;
t373 = -Ifges(5,3) - Ifges(6,2);
t400 = t373 * t322;
t205 = sin(pkin(8));
t206 = cos(pkin(8));
t209 = sin(qJ(2));
t282 = qJD(1) * t209;
t265 = t206 * t282;
t155 = qJD(2) * t205 + t265;
t208 = sin(qJ(4));
t211 = cos(qJ(4));
t263 = t205 * t282;
t230 = -t206 * qJD(2) + t263;
t219 = t211 * t155 - t208 * t230;
t331 = t219 / 0.2e1;
t375 = Ifges(6,4) + Ifges(5,5);
t399 = t375 * t331;
t194 = pkin(6) * t282;
t166 = -qJD(2) * pkin(2) + qJD(3) + t194;
t111 = pkin(3) * t230 + t166;
t255 = -qJ(3) * t209 - pkin(1);
t167 = -pkin(2) * t212 + t255;
t146 = t167 * qJD(1);
t195 = pkin(6) * t281;
t172 = qJD(2) * qJ(3) + t195;
t102 = t206 * t146 - t172 * t205;
t62 = -pkin(3) * t281 - pkin(7) * t155 + t102;
t103 = t205 * t146 + t206 * t172;
t65 = -pkin(7) * t230 + t103;
t19 = -t208 * t65 + t211 * t62;
t365 = -t19 + qJD(5);
t17 = -pkin(4) * t188 + t365;
t97 = t208 * t155 + t211 * t230;
t27 = t97 * pkin(4) - qJ(5) * t219 + t111;
t398 = mrSges(5,2) * t111 + t17 * mrSges(6,2) - t19 * mrSges(5,3) - mrSges(6,3) * t27;
t272 = qJD(1) * qJD(2);
t165 = qJDD(1) * t209 + t212 * t272;
t123 = qJDD(2) * t206 - t165 * t205;
t124 = qJDD(2) * t205 + t165 * t206;
t37 = -qJD(4) * t97 + t208 * t123 + t211 * t124;
t341 = t37 / 0.2e1;
t38 = qJD(4) * t219 - t211 * t123 + t208 * t124;
t339 = t38 / 0.2e1;
t333 = t97 / 0.2e1;
t384 = -m(5) - m(6);
t199 = t212 * qJDD(1);
t261 = t209 * t272;
t164 = -t199 + t261;
t157 = qJDD(4) + t164;
t325 = t157 / 0.2e1;
t258 = -t281 / 0.2e1;
t397 = -mrSges(5,3) - mrSges(6,2);
t377 = Ifges(5,1) + Ifges(6,1);
t374 = -Ifges(5,6) + Ifges(6,6);
t204 = pkin(8) + qJ(4);
t197 = sin(t204);
t198 = cos(t204);
t396 = -t378 * t197 + t198 * t379;
t316 = Ifges(6,5) * t97;
t95 = Ifges(5,4) * t97;
t372 = t375 * t188 + t377 * t219 + t316 - t95;
t395 = t372 / 0.2e1 + t398;
t210 = sin(qJ(1));
t213 = cos(qJ(1));
t393 = g(1) * t213 + g(2) * t210;
t20 = t208 * t62 + t211 * t65;
t18 = qJ(5) * t188 + t20;
t306 = Ifges(3,4) * t209;
t243 = t212 * Ifges(3,2) + t306;
t392 = t18 * mrSges(6,3) + t19 * mrSges(5,1) - Ifges(3,6) * qJD(2) / 0.2e1 - qJD(1) * t243 / 0.2e1 - t17 * mrSges(6,1) - t20 * mrSges(5,2) + Ifges(4,5) * t155 / 0.2e1 - Ifges(4,6) * t230 / 0.2e1 + Ifges(4,3) * t258 + t374 * t333 + t399 + t400;
t192 = pkin(6) * t199;
t127 = qJDD(2) * qJ(3) + t192 + (qJD(3) - t194) * qJD(2);
t278 = qJD(3) * t209;
t297 = qJDD(1) * pkin(1);
t88 = pkin(2) * t164 - qJ(3) * t165 - qJD(1) * t278 - t297;
t56 = -t127 * t205 + t206 * t88;
t26 = pkin(3) * t164 - pkin(7) * t124 + t56;
t57 = t206 * t127 + t205 * t88;
t29 = pkin(7) * t123 + t57;
t4 = -qJD(4) * t20 - t208 * t29 + t211 * t26;
t2 = -pkin(4) * t157 + qJDD(5) - t4;
t340 = -t38 / 0.2e1;
t154 = t165 * pkin(6);
t136 = -qJDD(2) * pkin(2) + qJDD(3) + t154;
t79 = -pkin(3) * t123 + t136;
t5 = pkin(4) * t38 - qJ(5) * t37 - qJD(5) * t219 + t79;
t391 = mrSges(5,2) * t79 + mrSges(6,2) * t2 - mrSges(5,3) * t4 - mrSges(6,3) * t5 + Ifges(5,4) * t340 + 0.2e1 * t325 * t375 + t339 * t401 + 0.2e1 * t341 * t377;
t94 = Ifges(6,5) * t219;
t41 = Ifges(6,6) * t188 + Ifges(6,3) * t97 + t94;
t317 = Ifges(5,4) * t219;
t44 = -Ifges(5,2) * t97 + Ifges(5,6) * t188 + t317;
t389 = t111 * mrSges(5,1) + t27 * mrSges(6,1) - mrSges(6,2) * t18 - mrSges(5,3) * t20 + t41 / 0.2e1 - t44 / 0.2e1;
t323 = t164 / 0.2e1;
t328 = t124 / 0.2e1;
t388 = Ifges(4,1) * t328 + Ifges(4,5) * t323;
t321 = t188 / 0.2e1;
t334 = -t97 / 0.2e1;
t387 = -Ifges(5,2) * t334 + Ifges(6,3) * t333 + t374 * t321 + t376 * t331 + t389;
t386 = Ifges(5,4) * t334 + Ifges(6,5) * t333 + t375 * t321 + t377 * t331 + t395;
t385 = -m(4) - m(3);
t329 = t123 / 0.2e1;
t382 = -t164 / 0.2e1;
t381 = t165 / 0.2e1;
t380 = qJD(2) / 0.2e1;
t159 = t205 * t211 + t206 * t208;
t228 = t159 * t212;
t121 = qJD(1) * t228;
t293 = t205 * t208;
t158 = -t211 * t206 + t293;
t227 = t158 * t212;
t122 = qJD(1) * t227;
t140 = t158 * qJD(4);
t141 = t159 * qJD(4);
t266 = t205 * t281;
t149 = pkin(3) * t266 + t195;
t371 = -qJD(5) * t159 - t149 + (-t122 + t140) * qJ(5) + (-t121 + t141) * pkin(4);
t70 = -mrSges(5,2) * t188 - mrSges(5,3) * t97;
t73 = -mrSges(6,2) * t97 + mrSges(6,3) * t188;
t370 = t70 + t73;
t318 = mrSges(5,3) * t219;
t71 = mrSges(5,1) * t188 - t318;
t319 = mrSges(6,2) * t219;
t72 = -mrSges(6,1) * t188 + t319;
t369 = -t72 + t71;
t250 = mrSges(3,1) * t212 - mrSges(3,2) * t209;
t368 = -mrSges(2,1) - t250;
t248 = -mrSges(4,1) * t206 + mrSges(4,2) * t205;
t226 = m(4) * pkin(2) - t248;
t367 = t212 * t226;
t366 = qJD(2) * mrSges(3,1) - mrSges(4,1) * t230 - t155 * mrSges(4,2) - mrSges(3,3) * t282;
t303 = Ifges(4,4) * t206;
t242 = -Ifges(4,2) * t205 + t303;
t304 = Ifges(4,4) * t205;
t244 = Ifges(4,1) * t206 - t304;
t364 = t155 * (Ifges(4,5) * t209 + t212 * t244) - (Ifges(4,6) * t209 + t212 * t242) * t230;
t363 = -t373 * t157 + t375 * t37 + t374 * t38;
t153 = -pkin(6) * t261 + t192;
t362 = t153 * t212 + t154 * t209;
t361 = t397 * t209;
t360 = -t205 * t56 + t206 * t57;
t247 = mrSges(4,1) * t205 + mrSges(4,2) * t206;
t358 = mrSges(2,2) - t247 - mrSges(3,3);
t289 = t206 * t212;
t291 = t205 * t212;
t357 = -t166 * t212 * t247 - t103 * (-mrSges(4,2) * t209 - mrSges(4,3) * t291) - t102 * (mrSges(4,1) * t209 - mrSges(4,3) * t289);
t356 = -t219 * t27 - t2;
t152 = t206 * t167;
t290 = t206 * t209;
t101 = -pkin(7) * t290 + t152 + (-pkin(6) * t205 - pkin(3)) * t212;
t118 = pkin(6) * t289 + t205 * t167;
t292 = t205 * t209;
t108 = -pkin(7) * t292 + t118;
t298 = t208 * t101 + t211 * t108;
t239 = pkin(2) * t209 - qJ(3) * t212;
t137 = qJD(2) * t239 - t278;
t280 = qJD(2) * t209;
t268 = pkin(6) * t280;
t106 = t206 * t137 + t205 * t268;
t233 = pkin(3) * t209 - pkin(7) * t289;
t74 = qJD(2) * t233 + t106;
t125 = t205 * t137;
t229 = -pkin(6) * t290 - pkin(7) * t291;
t83 = qJD(2) * t229 + t125;
t15 = -qJD(4) * t298 - t208 * t83 + t211 * t74;
t355 = m(6) * pkin(4) + t379;
t353 = m(6) * qJ(5) - t378;
t275 = qJD(4) * t211;
t277 = qJD(4) * t208;
t3 = t208 * t26 + t211 * t29 + t62 * t275 - t277 * t65;
t1 = qJ(5) * t157 + qJD(5) * t188 + t3;
t352 = t3 * mrSges(5,2) - t1 * mrSges(6,3);
t345 = mrSges(5,1) * t79 + mrSges(6,1) * t5 - mrSges(6,2) * t1 - mrSges(5,3) * t3 + 0.2e1 * Ifges(6,3) * t339 - t37 * Ifges(5,4) / 0.2e1 - t157 * Ifges(5,6) / 0.2e1 + t401 * t341 + (t374 + Ifges(6,6)) * t325 + (-t340 + t339) * Ifges(5,2);
t336 = Ifges(4,4) * t329 + t388;
t332 = -t219 / 0.2e1;
t315 = pkin(3) * t205;
t312 = g(3) * t209;
t200 = t209 * pkin(6);
t307 = pkin(7) + qJ(3);
t161 = t239 * qJD(1);
t142 = t205 * t161;
t100 = qJD(1) * t229 + t142;
t112 = pkin(6) * t263 + t206 * t161;
t82 = qJD(1) * t233 + t112;
t48 = t211 * t100 + t208 * t82;
t305 = Ifges(3,4) * t212;
t296 = t136 * t209;
t288 = t307 * t209;
t287 = t209 * t213;
t286 = t210 * t212;
t190 = pkin(3) * t206 + pkin(2);
t175 = t212 * t190;
t285 = t212 * t213;
t284 = t213 * t197;
t279 = qJD(2) * t212;
t196 = pkin(6) * t279;
t264 = t205 * t279;
t150 = pkin(3) * t264 + t196;
t162 = pkin(3) * t292 + t200;
t283 = t213 * pkin(1) + t210 * pkin(6);
t276 = qJD(4) * t209;
t262 = m(4) * qJ(3) + mrSges(4,3);
t11 = t38 * mrSges(5,1) + t37 * mrSges(5,2);
t10 = t38 * mrSges(6,1) - t37 * mrSges(6,3);
t23 = -t157 * mrSges(6,1) + t37 * mrSges(6,2);
t253 = -t272 / 0.2e1;
t252 = t272 / 0.2e1;
t67 = -t123 * mrSges(4,1) + t124 * mrSges(4,2);
t249 = mrSges(3,1) * t209 + mrSges(3,2) * t212;
t241 = Ifges(3,5) * t212 - Ifges(3,6) * t209;
t240 = Ifges(4,5) * t206 - Ifges(4,6) * t205;
t238 = pkin(4) * t198 + qJ(5) * t197;
t47 = -t100 * t208 + t211 * t82;
t54 = t101 * t211 - t108 * t208;
t169 = t307 * t205;
t170 = t307 * t206;
t234 = -t211 * t169 - t170 * t208;
t110 = -t169 * t208 + t170 * t211;
t232 = pkin(1) * t249;
t14 = t101 * t275 - t108 * t277 + t208 * t74 + t211 * t83;
t231 = t209 * (Ifges(3,1) * t212 - t306);
t220 = t212 * (Ifges(4,3) * t209 + t212 * t240);
t218 = t209 * t262 + t367;
t202 = t213 * pkin(6);
t193 = Ifges(3,4) * t281;
t173 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t281;
t145 = Ifges(3,1) * t282 + Ifges(3,5) * qJD(2) + t193;
t135 = t158 * t209;
t134 = t159 * t209;
t131 = t210 * t197 + t198 * t285;
t130 = -t210 * t198 + t212 * t284;
t129 = t198 * t286 - t284;
t128 = t197 * t286 + t198 * t213;
t120 = -mrSges(4,1) * t281 - mrSges(4,3) * t155;
t119 = mrSges(4,2) * t281 - mrSges(4,3) * t230;
t117 = -pkin(6) * t291 + t152;
t113 = -pkin(6) * t265 + t142;
t107 = -t206 * t268 + t125;
t93 = pkin(4) * t158 - qJ(5) * t159 - t190;
t91 = Ifges(4,1) * t155 - Ifges(4,4) * t230 - Ifges(4,5) * t281;
t90 = Ifges(4,4) * t155 - Ifges(4,2) * t230 - Ifges(4,6) * t281;
t86 = mrSges(4,1) * t164 - mrSges(4,3) * t124;
t85 = -mrSges(4,2) * t164 + mrSges(4,3) * t123;
t77 = qJD(2) * t228 + t275 * t290 - t276 * t293;
t76 = -qJD(2) * t227 - t159 * t276;
t69 = qJD(3) * t159 + qJD(4) * t110;
t68 = -qJD(3) * t158 + qJD(4) * t234;
t63 = pkin(4) * t134 + qJ(5) * t135 + t162;
t59 = t124 * Ifges(4,4) + t123 * Ifges(4,2) + t164 * Ifges(4,6);
t53 = mrSges(5,1) * t97 + mrSges(5,2) * t219;
t52 = mrSges(6,1) * t97 - mrSges(6,3) * t219;
t51 = pkin(4) * t219 + qJ(5) * t97;
t50 = pkin(4) * t212 - t54;
t49 = -qJ(5) * t212 + t298;
t40 = -pkin(4) * t282 - t47;
t39 = qJ(5) * t282 + t48;
t25 = -mrSges(6,2) * t38 + mrSges(6,3) * t157;
t24 = -mrSges(5,2) * t157 - mrSges(5,3) * t38;
t22 = mrSges(5,1) * t157 - mrSges(5,3) * t37;
t16 = pkin(4) * t77 - qJ(5) * t76 + qJD(5) * t135 + t150;
t13 = -pkin(4) * t280 - t15;
t12 = qJ(5) * t280 - qJD(5) * t212 + t14;
t6 = [-t391 * t135 + t387 * t77 + t386 * t76 + t250 * t297 + (Ifges(5,6) * t334 + Ifges(6,6) * t333 - t373 * t321 + t392 + t399) * t280 + t345 * t134 + (t206 * t91 + t145) * t279 / 0.2e1 + m(5) * (t111 * t150 + t14 * t20 + t15 * t19 + t162 * t79 + t298 * t3 + t4 * t54) + t298 * t24 - t366 * t196 + m(6) * (t1 * t49 + t12 * t18 + t13 * t17 + t16 * t27 + t2 * t50 + t5 * t63) + (-t290 * t56 - t292 * t57) * mrSges(4,3) + (t241 * t380 - t357) * qJD(2) + t290 * t336 + (-mrSges(3,1) * t200 + Ifges(3,5) * t209 + (-mrSges(3,2) * pkin(6) + Ifges(3,6)) * t212) * qJDD(2) + ((-Ifges(3,2) * t209 + t305) * t252 - Ifges(6,6) * t339 - Ifges(5,6) * t340 - Ifges(4,5) * t328 - Ifges(4,6) * t329 + Ifges(3,4) * t381 + Ifges(3,2) * t382 - Ifges(4,3) * t323 - t4 * mrSges(5,1) + t2 * mrSges(6,1) - t56 * mrSges(4,1) + t57 * mrSges(4,2) - t375 * t341 + t373 * t325 + t352) * t212 + (Ifges(3,1) * t165 + Ifges(3,4) * t382 + t240 * t323 + t242 * t329 + t244 * t328) * t209 + (t384 * (-t210 * t288 + t213 * t315 + t202) + t385 * t202 + t355 * t129 + t353 * t128 + t358 * t213 + (m(3) * pkin(1) - m(4) * t255 + t209 * mrSges(4,3) + t367 + t384 * (-pkin(1) - t175) - t361 - t368) * t210) * g(1) + t231 * t252 + t220 * t253 + t364 * t380 + t305 * t381 + t243 * t382 + t49 * t25 + t50 * t23 + t16 * t52 + t54 * t22 + t63 * t10 + t14 * t70 + t15 * t71 + t13 * t72 + t12 * t73 + t117 * t86 + t118 * t85 + t107 * t119 + t106 * t120 + t247 * t296 + t67 * t200 + t150 * t53 + t162 * t11 - pkin(1) * (mrSges(3,1) * t164 + mrSges(3,2) * t165) + (t397 * t287 + t385 * t283 + t384 * (t190 * t285 + t210 * t315 + t287 * t307 + t283) - t355 * t131 - t353 * t130 + (-t218 + t368) * t213 + t358 * t210) * g(2) + (-pkin(6) * t164 * t212 + t165 * t200 + t362) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(6) * t362) - (Ifges(4,5) * t124 + Ifges(4,6) * t123 + Ifges(4,3) * t164 + t363) * t212 / 0.2e1 - t90 * t264 / 0.2e1 - t173 * t268 + Ifges(2,3) * qJDD(1) - t232 * t272 - t59 * t292 / 0.2e1 + m(4) * (t102 * t106 + t103 * t107 + t117 * t56 + t118 * t57 + (t166 * t279 + t296) * pkin(6)); t391 * t159 + (-Ifges(5,2) * t333 + Ifges(6,3) * t334 + t374 * t322 + t376 * t332 - t389) * t121 + (-qJ(3) * t86 + t336 + (-m(4) * t102 - t120) * qJD(3) + t388) * t205 + t387 * t141 - t386 * t140 + (t234 * t4 + t110 * t3 - t111 * t149 - t190 * t79 + (-t48 + t68) * t20 + (-t47 - t69) * t19) * m(5) + (t1 * t110 - t234 * t2 + t5 * t93 + t371 * t27 + (-t39 + t68) * t18 + (-t40 + t69) * t17) * m(6) - (t23 - t22) * t234 + (-Ifges(3,2) * t258 + Ifges(5,6) * t333 + Ifges(6,6) * t334 + t375 * t332 - t392 - t400) * t282 + t304 * t329 + t345 * t158 + t366 * t195 + (t145 + t193) * t258 + (t357 - t364 / 0.2e1 + (t220 / 0.2e1 - t231 / 0.2e1 + t232) * qJD(1)) * qJD(1) - (Ifges(5,4) * t333 + Ifges(6,5) * t334 + t375 * t322 + t377 * t332 - t395) * t122 + (-t250 - t218 + t384 * (t175 + t288) + (-m(6) * t238 - t396) * t212 + t361) * g(3) + t303 * t328 - t369 * t69 + t370 * t68 + t371 * t52 + (t91 * t258 + Ifges(4,2) * t329 + Ifges(4,6) * t323 + qJ(3) * t85 + t59 / 0.2e1 + (m(4) * t103 + t119) * qJD(3)) * t206 + t173 * t194 + t241 * t253 + t393 * (t249 + (t384 * t307 - t262 + t397) * t212 + (t226 + m(5) * t190 - m(6) * (-t190 - t238) + t396) * t209) - pkin(2) * t67 - t48 * t70 - t47 * t71 - t40 * t72 - t39 * t73 + t136 * t248 + t93 * t10 - t113 * t119 - t112 * t120 - t149 * t53 - t153 * mrSges(3,2) - t154 * mrSges(3,1) - Ifges(3,6) * t164 + Ifges(3,5) * t165 + Ifges(3,3) * qJDD(2) - t190 * t11 + t360 * mrSges(4,3) + (-pkin(2) * t136 + qJ(3) * t360 - t102 * t112 - t103 * t113 - t166 * t195) * m(4) + t90 * t266 / 0.2e1 + (t24 + t25) * t110; t155 * t120 + t370 * t97 + t369 * t219 + t67 + t10 + t11 + t230 * t119 + (t212 * g(3) - t209 * t393) * (m(4) - t384) + (-t17 * t219 + t18 * t97 + t5) * m(6) + (t19 * t219 + t20 * t97 + t79) * m(5) + (t102 * t155 + t103 * t230 + t136) * m(4); (-m(6) * t17 + t318 + t369) * t20 + (Ifges(6,3) * t219 - t316) * t334 + (-t111 * t219 + t4) * mrSges(5,1) + (t219 * t374 - t375 * t97) * t322 + (-Ifges(5,2) * t219 + t372 - t95) * t333 - t352 + (t379 * t197 + t378 * t198) * t312 + (t379 * t128 + t378 * t129) * g(2) + (t130 * t379 + t131 * t378) * g(1) + t44 * t331 + (-t377 * t97 - t317 + t41 + t94) * t332 + (-t27 * t51 - (-pkin(4) * t197 + qJ(5) * t198) * t312 - pkin(4) * t2 + qJ(5) * t1 - g(2) * (-pkin(4) * t128 + qJ(5) * t129) - g(1) * (-pkin(4) * t130 + qJ(5) * t131) + t365 * t18) * m(6) + t356 * mrSges(6,1) + t398 * t97 + t18 * t319 - pkin(4) * t23 + qJ(5) * t25 - t51 * t52 + t363 + qJD(5) * t73 - t370 * t19; -t188 * t73 + t219 * t52 + (-g(1) * t130 - g(2) * t128 - t18 * t188 - t197 * t312 - t356) * m(6) + t23;];
tau = t6;
