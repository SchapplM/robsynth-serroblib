% Calculate vector of inverse dynamics joint torques for
% S6PRPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% tau [6x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PRPRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:24:45
% EndTime: 2019-03-08 19:25:07
% DurationCPUTime: 14.25s
% Computational Cost: add. (5935->592), mult. (13602->828), div. (0->0), fcn. (10750->16), ass. (0->271)
t309 = cos(pkin(12));
t193 = sin(pkin(12));
t205 = cos(qJ(4));
t196 = sin(pkin(6));
t206 = cos(qJ(2));
t279 = qJD(1) * t206;
t255 = t196 * t279;
t162 = qJD(2) * pkin(2) + t255;
t194 = sin(pkin(11));
t197 = cos(pkin(11));
t203 = sin(qJ(2));
t290 = t196 * t203;
t256 = qJD(1) * t290;
t116 = t162 * t194 + t197 * t256;
t112 = qJD(2) * pkin(8) + t116;
t246 = qJ(5) * qJD(2) + t112;
t199 = cos(pkin(6));
t179 = qJD(1) * t199 + qJD(3);
t202 = sin(qJ(4));
t296 = t179 * t202;
t71 = t205 * t246 + t296;
t319 = t193 * t71;
t166 = t205 * t179;
t70 = -t202 * t246 + t166;
t67 = qJD(4) * pkin(4) + t70;
t27 = t309 * t67 - t319;
t250 = t309 * t202;
t277 = qJD(2) * t205;
t146 = -qJD(2) * t250 - t193 * t277;
t201 = sin(qJ(6));
t204 = cos(qJ(6));
t119 = qJD(4) * t204 + t146 * t201;
t120 = qJD(4) * t201 - t146 * t204;
t326 = mrSges(6,3) * t146;
t310 = -qJD(4) * mrSges(6,1) - mrSges(7,1) * t119 + mrSges(7,2) * t120 - t326;
t381 = m(6) * t27 - t310;
t237 = -mrSges(7,1) * t204 + mrSges(7,2) * t201;
t379 = m(7) * pkin(5) + mrSges(6,1) - t237;
t185 = pkin(2) * t194 + pkin(8);
t281 = qJ(5) + t185;
t247 = qJD(4) * t281;
t126 = qJD(5) * t205 - t202 * t247;
t165 = t194 * t256;
t130 = t197 * t255 - t165;
t211 = -qJD(5) * t202 - t205 * t247;
t217 = -t193 * t202 + t205 * t309;
t378 = -t126 * t309 + t130 * t217 - t193 * t211;
t228 = t194 * t206 + t197 * t203;
t137 = t228 * t196;
t127 = qJD(1) * t137;
t155 = t193 * t205 + t250;
t147 = t155 * qJD(4);
t148 = t217 * qJD(4);
t276 = qJD(4) * t202;
t268 = pkin(4) * t276;
t377 = pkin(5) * t147 - pkin(9) * t148 - t127 + t268;
t371 = -m(6) - m(7);
t372 = -m(5) - m(4);
t376 = t371 + t372;
t192 = qJ(4) + pkin(12);
t190 = sin(t192);
t191 = cos(t192);
t239 = -mrSges(5,1) * t205 + mrSges(5,2) * t202;
t375 = m(5) * pkin(3) - t239 + t379 * t191 + (m(7) * pkin(9) - mrSges(6,2) + mrSges(7,3)) * t190;
t195 = sin(pkin(10));
t198 = cos(pkin(10));
t283 = t199 * t206;
t374 = -t195 * t203 + t198 * t283;
t289 = t196 * t205;
t282 = t206 * t197;
t154 = t194 * t203 - t282;
t285 = t199 * t203;
t280 = -t194 * t283 - t197 * t285;
t91 = -t154 * t198 + t195 * t280;
t373 = t195 * t289 - t202 * t91;
t109 = -t137 * t202 + t199 * t205;
t145 = t217 * qJD(2);
t142 = qJD(6) - t145;
t272 = qJD(2) * qJD(4);
t159 = qJDD(2) * t205 - t202 * t272;
t160 = qJDD(2) * t202 + t205 * t272;
t108 = t159 * t193 + t160 * t309;
t47 = qJD(6) * t119 + qJDD(4) * t201 + t108 * t204;
t344 = t47 / 0.2e1;
t48 = -qJD(6) * t120 + qJDD(4) * t204 - t108 * t201;
t343 = t48 / 0.2e1;
t64 = t309 * t71;
t28 = t193 * t67 + t64;
t26 = qJD(4) * pkin(9) + t28;
t115 = t162 * t197 - t165;
t188 = pkin(4) * t205 + pkin(3);
t99 = -qJD(2) * t188 + qJD(5) - t115;
t44 = -pkin(5) * t145 + pkin(9) * t146 + t99;
t9 = -t201 * t26 + t204 * t44;
t370 = t9 * mrSges(7,1);
t107 = t159 * t309 - t160 * t193;
t106 = qJDD(6) - t107;
t342 = t106 / 0.2e1;
t369 = -t147 / 0.2e1;
t368 = t148 / 0.2e1;
t10 = t201 * t44 + t204 * t26;
t367 = t10 * mrSges(7,2);
t333 = pkin(2) * t197;
t167 = -t188 - t333;
t92 = -pkin(5) * t217 - pkin(9) * t155 + t167;
t152 = t281 * t205;
t248 = t281 * t202;
t96 = t152 * t309 - t193 * t248;
t39 = -t201 * t96 + t204 * t92;
t366 = qJD(6) * t39 + t201 * t377 - t204 * t378;
t40 = t201 * t92 + t204 * t96;
t365 = -qJD(6) * t40 + t201 * t378 + t204 * t377;
t17 = -mrSges(7,1) * t48 + mrSges(7,2) * t47;
t98 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t108;
t328 = t17 - t98;
t311 = qJDD(4) / 0.2e1;
t136 = t194 * t290 - t196 * t282;
t125 = t136 * t204;
t110 = t137 * t205 + t199 * t202;
t50 = t109 * t193 + t110 * t309;
t35 = -t201 * t50 + t125;
t360 = (mrSges(3,1) * t206 - mrSges(3,2) * t203) * t196 - t136 * mrSges(4,1) - t137 * mrSges(4,2);
t273 = qJD(6) * t204;
t222 = t148 * t201 + t155 * t273;
t86 = t154 * t195 + t198 * t280;
t359 = -t198 * t289 + t202 * t86;
t278 = qJD(2) * t202;
t265 = mrSges(5,3) * t278;
t170 = qJD(4) * mrSges(5,1) - t265;
t264 = mrSges(5,3) * t277;
t171 = -qJD(4) * mrSges(5,2) + t264;
t358 = -t170 * t202 + t171 * t205;
t178 = qJDD(1) * t199 + qJDD(3);
t275 = qJD(4) * t205;
t244 = qJD(2) * t256;
t288 = t196 * t206;
t140 = qJDD(1) * t288 - t244;
t306 = qJDD(2) * pkin(2);
t135 = t140 + t306;
t253 = qJD(2) * t279;
t141 = (qJDD(1) * t203 + t253) * t196;
t77 = t135 * t194 + t141 * t197;
t73 = qJDD(2) * pkin(8) + t77;
t29 = -t112 * t276 + t178 * t202 + t179 * t275 + t205 * t73;
t164 = t205 * t178;
t83 = t112 * t205 + t296;
t30 = -qJD(4) * t83 - t202 * t73 + t164;
t357 = -t202 * t30 + t205 * t29;
t33 = mrSges(7,1) * t106 - mrSges(7,3) * t47;
t34 = -mrSges(7,2) * t106 + mrSges(7,3) * t48;
t356 = -t201 * t33 + t204 * t34;
t76 = t135 * t197 - t141 * t194;
t72 = -qJDD(2) * pkin(3) - t76;
t52 = -pkin(4) * t159 + qJDD(5) + t72;
t21 = -pkin(5) * t107 - pkin(9) * t108 + t52;
t271 = qJD(2) * qJD(5);
t20 = -t112 * t275 + qJDD(4) * pkin(4) - qJ(5) * t160 + t164 + (-qJD(4) * t179 - t271 - t73) * t202;
t22 = qJ(5) * t159 + t205 * t271 + t29;
t6 = t193 * t20 + t22 * t309;
t4 = qJDD(4) * pkin(9) + t6;
t1 = qJD(6) * t9 + t201 * t21 + t204 * t4;
t2 = -qJD(6) * t10 - t201 * t4 + t204 * t21;
t355 = t1 * t204 - t2 * t201;
t354 = 0.2e1 * t311;
t100 = -mrSges(6,1) * t145 - mrSges(6,2) * t146;
t353 = t100 + (-mrSges(4,1) + t239) * qJD(2);
t352 = -mrSges(4,1) - t375;
t351 = t2 * mrSges(7,1) - t1 * mrSges(7,2);
t236 = mrSges(7,1) * t201 + mrSges(7,2) * t204;
t350 = -m(5) * pkin(8) - mrSges(5,3) - mrSges(6,3) - t236;
t349 = -mrSges(4,2) - t350;
t348 = qJD(2) ^ 2;
t347 = m(6) * pkin(4);
t346 = Ifges(7,1) * t344 + Ifges(7,4) * t343 + Ifges(7,5) * t342;
t341 = -t119 / 0.2e1;
t340 = -t120 / 0.2e1;
t339 = t120 / 0.2e1;
t338 = -t142 / 0.2e1;
t336 = -t146 / 0.2e1;
t334 = t204 / 0.2e1;
t332 = pkin(4) * t193;
t327 = mrSges(6,3) * t145;
t325 = Ifges(5,4) * t202;
t324 = Ifges(5,4) * t205;
t323 = Ifges(6,4) * t146;
t322 = Ifges(7,4) * t201;
t321 = Ifges(7,4) * t204;
t320 = t120 * Ifges(7,4);
t308 = mrSges(4,2) * qJD(2);
t305 = t136 * t201;
t303 = t145 * t201;
t302 = t145 * t204;
t300 = t155 * t201;
t299 = t155 * t204;
t294 = t195 * t196;
t292 = t196 * t198;
t291 = t196 * t202;
t274 = qJD(6) * t201;
t180 = pkin(2) * t288;
t270 = Ifges(7,5) * t47 + Ifges(7,6) * t48 + Ifges(7,3) * t106;
t269 = pkin(4) * t278;
t114 = Ifges(7,4) * t119;
t43 = Ifges(7,1) * t120 + Ifges(7,5) * t142 + t114;
t259 = t43 * t334;
t257 = t309 * pkin(4);
t251 = -t274 / 0.2e1;
t51 = -mrSges(6,1) * t107 + mrSges(6,2) * t108;
t245 = t373 * pkin(4);
t243 = t374 * pkin(2);
t242 = t109 * pkin(4);
t241 = t10 * t204 - t201 * t9;
t235 = Ifges(7,1) * t204 - t322;
t234 = t205 * Ifges(5,2) + t325;
t233 = -Ifges(7,2) * t201 + t321;
t232 = Ifges(5,5) * t205 - Ifges(5,6) * t202;
t231 = Ifges(7,5) * t204 - Ifges(7,6) * t201;
t78 = -mrSges(7,2) * t142 + mrSges(7,3) * t119;
t79 = mrSges(7,1) * t142 - mrSges(7,3) * t120;
t230 = -t201 * t79 + t204 * t78;
t82 = -t112 * t202 + t166;
t229 = -t202 * t82 + t205 * t83;
t36 = t204 * t50 + t305;
t131 = -qJD(4) * mrSges(6,2) + t327;
t227 = -t131 - t230;
t25 = -qJD(4) * pkin(5) - t27;
t223 = t25 * t236;
t5 = -t193 * t22 + t20 * t309;
t221 = -t148 * t204 + t155 * t274;
t111 = -qJD(2) * pkin(3) - t115;
t220 = t111 * (mrSges(5,1) * t202 + mrSges(5,2) * t205);
t219 = t202 * (Ifges(5,1) * t205 - t325);
t215 = t154 * t199;
t87 = -t195 * t228 - t198 * t215;
t90 = t195 * t215 - t198 * t228;
t216 = g(1) * t90 + g(2) * t87 - g(3) * t136;
t209 = (-t10 * t201 - t204 * t9) * qJD(6) + t355;
t200 = -qJ(5) - pkin(8);
t189 = Ifges(5,4) * t277;
t186 = -t257 - pkin(5);
t150 = Ifges(5,1) * t278 + Ifges(5,5) * qJD(4) + t189;
t149 = Ifges(5,6) * qJD(4) + qJD(2) * t234;
t144 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t160;
t143 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t159;
t139 = Ifges(6,4) * t145;
t129 = t154 * t196 * qJD(2);
t128 = qJD(2) * t137;
t113 = -mrSges(5,1) * t159 + mrSges(5,2) * t160;
t103 = t137 * t191 + t190 * t199;
t97 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t107;
t95 = t152 * t193 + t248 * t309;
t94 = -t146 * Ifges(6,1) + Ifges(6,5) * qJD(4) + t139;
t93 = t145 * Ifges(6,2) + Ifges(6,6) * qJD(4) - t323;
t84 = -pkin(5) * t146 - pkin(9) * t145 + t269;
t59 = qJD(4) * t109 - t129 * t205;
t58 = -qJD(4) * t110 + t129 * t202;
t57 = t190 * t294 + t191 * t91;
t55 = -t190 * t292 - t191 * t86;
t49 = -t109 * t309 + t110 * t193;
t42 = t119 * Ifges(7,2) + t142 * Ifges(7,6) + t320;
t41 = Ifges(7,5) * t120 + t119 * Ifges(7,6) + t142 * Ifges(7,3);
t32 = t309 * t70 - t319;
t31 = t193 * t70 + t64;
t24 = t193 * t58 + t309 * t59;
t23 = t193 * t59 - t309 * t58;
t16 = t201 * t84 + t204 * t32;
t15 = -t201 * t32 + t204 * t84;
t13 = Ifges(7,4) * t47 + Ifges(7,2) * t48 + Ifges(7,6) * t106;
t8 = -qJD(6) * t36 + t128 * t204 - t201 * t24;
t7 = qJD(6) * t35 + t128 * t201 + t204 * t24;
t3 = -qJDD(4) * pkin(5) - t5;
t11 = [t129 * t308 + m(2) * qJDD(1) + t109 * t144 + t110 * t143 + t24 * t131 + t58 * t170 + t59 * t171 + t35 * t33 + t36 * t34 + t50 * t97 + t7 * t78 + t8 * t79 + t328 * t49 + t310 * t23 + (-mrSges(3,1) * t203 - mrSges(3,2) * t206) * t348 * t196 + (t113 + t51) * t136 + t353 * t128 + t360 * qJDD(2) + (-m(2) - m(3) + t376) * g(3) + m(7) * (t1 * t36 + t10 * t7 + t2 * t35 + t23 * t25 + t3 * t49 + t8 * t9) + m(6) * (t128 * t99 + t136 * t52 - t23 * t27 + t24 * t28 - t49 * t5 + t50 * t6) + m(5) * (t109 * t30 + t110 * t29 + t111 * t128 + t136 * t72 + t58 * t82 + t59 * t83) + m(4) * (-t115 * t128 - t116 * t129 - t136 * t76 + t137 * t77 + t178 * t199) + m(3) * (qJDD(1) * t199 ^ 2 + (t140 * t206 + t141 * t203) * t196); (-t374 * mrSges(3,1) - (-t195 * t206 - t198 * t285) * mrSges(3,2) + t372 * t243 + t371 * (t188 * t87 + t200 * t86 + t243) + t352 * t87 + t349 * t86) * g(2) + (t371 * (-t136 * t188 - t137 * t200 + t180) + t372 * t180 + t350 * t137 + t375 * t136 - t360) * g(3) + (-m(7) * t25 + t381) * (-t126 * t193 + t130 * t155 + t211 * t309) + t159 * t234 / 0.2e1 + (t1 * t40 + t10 * t366 + t2 * t39 + t3 * t95 + t365 * t9) * m(7) + (t220 + Ifges(6,5) * t368 + Ifges(6,6) * t369 + t232 * qJD(4) / 0.2e1) * qJD(4) + (m(5) * ((-t202 * t83 - t205 * t82) * qJD(4) + t357) - t170 * t275 - t171 * t276 - t202 * t144 + t205 * t143) * t185 + (-t275 * t82 - t276 * t83 + t357) * mrSges(5,3) + (-m(4) * t116 - m(5) * t229 + t308 - t358) * t130 - t222 * t42 / 0.2e1 + t25 * (mrSges(7,1) * t222 - mrSges(7,2) * t221) + t142 * (-Ifges(7,5) * t221 - Ifges(7,6) * t222 + Ifges(7,3) * t147) / 0.2e1 + t119 * (-Ifges(7,4) * t221 - Ifges(7,2) * t222 + Ifges(7,6) * t147) / 0.2e1 + t299 * t346 + (Ifges(6,1) * t148 - Ifges(6,4) * t147) * t336 + (-Ifges(7,1) * t221 - Ifges(7,4) * t222 + Ifges(7,5) * t147) * t339 + t94 * t368 + t93 * t369 + t147 * t370 + (Ifges(5,5) * t202 + Ifges(5,6) * t205) * t311 + t148 * t259 - t147 * t367 + t160 * t202 * Ifges(5,1) + (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (-t1 * t300 - t10 * t222 - t2 * t299 + t221 * t9) * mrSges(7,3) + t328 * t95 + t365 * t79 + t366 * t78 - (Ifges(7,3) * t342 + Ifges(7,6) * t343 + Ifges(7,5) * t344 + t270 / 0.2e1 - Ifges(6,2) * t107 - Ifges(6,4) * t108 + t52 * mrSges(6,1) - t6 * mrSges(6,3) - t354 * Ifges(6,6) + t351) * t217 + m(4) * (t194 * t77 + t197 * t76) * pkin(2) + (-(t195 * t285 - t198 * t206) * mrSges(3,2) + t371 * (t188 * t90 - t200 * t91) + t352 * t90 - t349 * t91 + (pkin(2) * t376 - mrSges(3,1)) * (-t195 * t283 - t198 * t203)) * g(1) - t13 * t300 / 0.2e1 + (t167 * t52 + t268 * t99 - t28 * t378 - t5 * t95 + t6 * t96) * m(6) - t378 * t131 + (m(4) * t115 - m(5) * t111 - m(6) * t99 - t353) * t127 + (t52 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,1) * t108 + Ifges(6,4) * t107 + Ifges(6,5) * t354 + t231 * t342 + t233 * t343 + t235 * t344 + t3 * t236 + t43 * t251) * t155 + t160 * t324 / 0.2e1 + t202 * (Ifges(5,4) * t159 + Ifges(5,5) * qJDD(4)) / 0.2e1 + (-t194 * t306 - t77) * mrSges(4,2) + (t196 * t253 - t141) * mrSges(3,2) + (t197 * t306 + t76) * mrSges(4,1) + (-t147 * t28 - t148 * t27) * mrSges(6,3) + (t244 + t140) * mrSges(3,1) + t205 * (Ifges(5,4) * t160 + Ifges(5,2) * t159 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t167 * t51 + t147 * t41 / 0.2e1 + t145 * (Ifges(6,4) * t148 - Ifges(6,2) * t147) / 0.2e1 + t99 * (mrSges(6,1) * t147 + mrSges(6,2) * t148) + t96 * t97 + t150 * t275 / 0.2e1 - t149 * t276 / 0.2e1 + t100 * t268 + t39 * t33 + t40 * t34 + t72 * t239 + (m(5) * t72 + t113) * (-pkin(3) - t333) + (t219 + t205 * (-Ifges(5,2) * t202 + t324)) * t272 / 0.2e1; t202 * t143 + t205 * t144 - t328 * t217 + t310 * t147 + t358 * qJD(4) - t227 * t148 + (t97 + (-t201 * t78 - t204 * t79) * qJD(6) + t356) * t155 + m(6) * (-t147 * t27 + t148 * t28 + t155 * t6 + t217 * t5) + m(4) * t178 + m(7) * (t147 * t25 + t148 * t241 + t155 * t209 - t217 * t3) + m(5) * (qJD(4) * t229 + t202 * t29 + t205 * t30) - (-t199 * g(3) + (-g(1) * t195 + g(2) * t198) * t196) * t376; t381 * t31 + t3 * t237 - t145 * t223 - m(6) * (t269 * t99 + t28 * t32) + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) - t43 * t302 / 0.2e1 - qJD(2) * t220 + (-Ifges(7,5) * t146 + t145 * t235) * t340 + (-Ifges(7,6) * t146 + t145 * t233) * t341 + (Ifges(7,5) * t201 + Ifges(7,6) * t204) * t342 + (Ifges(7,2) * t204 + t322) * t343 + (Ifges(7,1) * t201 + t321) * t344 + t201 * t346 + (t193 * t6 + t309 * t5) * t347 + t27 * t327 + t97 * t332 + t13 * t334 + t93 * t336 + (-Ifges(7,3) * t146 + t145 * t231) * t338 + t146 * t370 + t98 * t257 + (t303 / 0.2e1 + t251) * t42 - t146 * t367 + t149 * t278 / 0.2e1 - t100 * t269 + (-t10 * t16 - t15 * t9 + t186 * t3 - t25 * t31) * m(7) - t348 * t219 / 0.2e1 - t28 * t326 + (-g(1) * t57 - g(2) * t55 - g(3) * t103 + (-t273 + t302) * t9 + (-t274 + t303) * t10 + t355) * mrSges(7,3) + (m(7) * t209 - t273 * t79 - t274 * t78 + t356) * (pkin(9) + t332) + (-t171 + t264) * t82 + (t223 + t259) * qJD(6) + (-m(6) * t245 + t57 * mrSges(6,2) - t373 * mrSges(5,1) - (-t195 * t291 - t205 * t91) * mrSges(5,2) - m(7) * (pkin(9) * t57 + t245) - t379 * (-t190 * t91 + t191 * t294)) * g(1) + (-m(6) * t242 + t103 * mrSges(6,2) - mrSges(5,1) * t109 + mrSges(5,2) * t110 - m(7) * (pkin(9) * t103 + t242) - t379 * (-t137 * t190 + t191 * t199)) * g(3) + (-(t198 * t291 + t205 * t86) * mrSges(5,2) + t55 * mrSges(6,2) - m(7) * (pkin(4) * t359 + pkin(9) * t55) - t379 * (t190 * t86 - t191 * t292) + (-mrSges(5,1) - t347) * t359) * g(2) + t186 * t17 + Ifges(5,6) * t159 + Ifges(5,5) * t160 - qJD(4) * (Ifges(6,5) * t145 + Ifges(6,6) * t146) / 0.2e1 - t99 * (-mrSges(6,1) * t146 + mrSges(6,2) * t145) - t32 * t131 + Ifges(6,6) * t107 + Ifges(6,5) * t108 + (t170 + t265) * t83 - t232 * t272 / 0.2e1 + t5 * mrSges(6,1) - t6 * mrSges(6,2) - t29 * mrSges(5,2) + t30 * mrSges(5,1) - t16 * t78 - t15 * t79 - (-Ifges(5,2) * t278 + t150 + t189) * t277 / 0.2e1 + (t119 * t233 + t120 * t235 + t142 * t231) * qJD(6) / 0.2e1 + (Ifges(6,1) * t145 + t323 + t41) * t146 / 0.2e1 - (Ifges(6,2) * t146 + t139 + t94) * t145 / 0.2e1; t230 * qJD(6) + t227 * t145 + t310 * t146 + t201 * t34 + t204 * t33 + t51 + (t1 * t201 + t142 * t241 + t146 * t25 + t2 * t204 + t216) * m(7) + (-t145 * t28 - t146 * t27 + t216 + t52) * m(6); -t25 * (mrSges(7,1) * t120 + mrSges(7,2) * t119) + (Ifges(7,1) * t119 - t320) * t340 + t42 * t339 + (Ifges(7,5) * t119 - Ifges(7,6) * t120) * t338 + t10 * t79 - t9 * t78 - g(1) * ((-t201 * t57 - t204 * t90) * mrSges(7,1) + (t201 * t90 - t204 * t57) * mrSges(7,2)) - g(2) * ((-t201 * t55 - t204 * t87) * mrSges(7,1) + (t201 * t87 - t204 * t55) * mrSges(7,2)) - g(3) * ((-t103 * t201 + t125) * mrSges(7,1) + (-t103 * t204 - t305) * mrSges(7,2)) + (t10 * t120 + t119 * t9) * mrSges(7,3) + t270 + (-Ifges(7,2) * t120 + t114 + t43) * t341 + t351;];
tau  = t11;
