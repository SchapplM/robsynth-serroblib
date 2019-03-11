% Calculate vector of inverse dynamics joint torques for
% S6PPRPRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
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
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6PPRPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:42:04
% EndTime: 2019-03-08 18:42:22
% DurationCPUTime: 9.20s
% Computational Cost: add. (5907->506), mult. (15448->749), div. (0->0), fcn. (14235->16), ass. (0->242)
t176 = cos(pkin(6));
t170 = sin(pkin(7));
t182 = cos(qJ(3));
t269 = t170 * t182;
t171 = sin(pkin(6));
t168 = sin(pkin(12));
t173 = cos(pkin(12));
t179 = sin(qJ(3));
t175 = cos(pkin(7));
t265 = t175 * t182;
t204 = -t168 * t179 + t173 * t265;
t328 = t171 * t204;
t97 = t176 * t269 + t328;
t334 = m(6) + m(7);
t178 = sin(qJ(5));
t181 = cos(qJ(5));
t224 = pkin(5) * t178 - pkin(10) * t181;
t167 = sin(pkin(13));
t172 = cos(pkin(13));
t159 = qJD(1) * t176 + qJD(2);
t270 = t170 * t179;
t203 = t173 * t175 * t179 + t168 * t182;
t320 = t203 * t171;
t86 = qJD(1) * t320 + t159 * t270;
t81 = t172 * t86;
t85 = qJD(1) * t328 + t159 * t269;
t55 = t167 * t85 + t81;
t333 = t224 * qJD(5) - t55;
t177 = sin(qJ(6));
t180 = cos(qJ(6));
t254 = qJD(5) * t180;
t258 = qJD(3) * t178;
t146 = -t177 * t258 + t254;
t332 = -t146 / 0.2e1;
t256 = qJD(5) * t177;
t147 = t180 * t258 + t256;
t331 = -t147 / 0.2e1;
t257 = qJD(3) * t181;
t160 = qJD(6) - t257;
t330 = -t160 / 0.2e1;
t244 = mrSges(6,3) * t258;
t259 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t146 + mrSges(7,2) * t147 + t244;
t271 = t170 * t171;
t239 = t173 * t271;
t114 = -qJD(1) * t239 + t159 * t175 + qJD(4);
t83 = qJD(3) * pkin(3) + t85;
t53 = t167 * t83 + t81;
t49 = qJD(3) * pkin(9) + t53;
t34 = t114 * t178 + t181 * t49;
t32 = qJD(5) * pkin(10) + t34;
t207 = -pkin(5) * t181 - pkin(10) * t178 - pkin(4);
t80 = t167 * t86;
t52 = t172 * t83 - t80;
t35 = qJD(3) * t207 - t52;
t10 = t177 * t35 + t180 * t32;
t249 = qJD(3) * qJD(5);
t151 = qJDD(3) * t181 - t178 * t249;
t152 = qJDD(3) * t178 + t181 * t249;
t280 = qJDD(3) * pkin(3);
t158 = qJDD(1) * t176 + qJDD(2);
t61 = -t86 * qJD(3) + qJDD(1) * t328 + t158 * t269;
t59 = t61 + t280;
t60 = (qJD(3) * t159 * t182 + t158 * t179) * t170 + (qJD(1) * qJD(3) * t204 + qJDD(1) * t203) * t171;
t21 = -t167 * t60 + t172 * t59;
t17 = -qJDD(3) * pkin(4) - t21;
t11 = -pkin(5) * t151 - pkin(10) * t152 + t17;
t117 = -qJDD(1) * t239 + t175 * t158;
t111 = qJDD(4) + t117;
t22 = t167 * t59 + t172 * t60;
t18 = qJDD(3) * pkin(9) + t22;
t253 = qJD(5) * t181;
t255 = qJD(5) * t178;
t7 = t178 * t111 + t114 * t253 + t181 * t18 - t255 * t49;
t5 = qJDD(5) * pkin(10) + t7;
t9 = -t177 * t32 + t180 * t35;
t1 = qJD(6) * t9 + t11 * t177 + t180 * t5;
t2 = -qJD(6) * t10 + t11 * t180 - t177 * t5;
t223 = t1 * t180 - t177 * t2;
t250 = qJD(6) * t180;
t252 = qJD(6) * t177;
t329 = -t10 * t252 - t9 * t250 + t223;
t98 = t176 * t270 + t320;
t326 = t151 / 0.2e1;
t325 = t152 / 0.2e1;
t300 = pkin(3) * t172;
t135 = t207 - t300;
t163 = pkin(3) * t167 + pkin(9);
t261 = t180 * t181;
t102 = t135 * t177 + t163 * t261;
t234 = t163 * t255;
t263 = t177 * t181;
t56 = t172 * t85 - t80;
t324 = -qJD(6) * t102 + t177 * t234 + t180 * t333 + t263 * t56;
t101 = t135 * t180 - t163 * t263;
t323 = qJD(6) * t101 + t177 * t333 - t180 * t234 - t261 * t56;
t221 = -mrSges(7,1) * t180 + mrSges(7,2) * t177;
t196 = m(7) * pkin(5) - t221;
t322 = -mrSges(6,1) - t196;
t245 = m(7) * pkin(10) + mrSges(7,3);
t321 = mrSges(6,2) - t245;
t94 = qJD(6) * t146 + qJDD(5) * t177 + t152 * t180;
t95 = -qJD(6) * t147 + qJDD(5) * t180 - t152 * t177;
t64 = -mrSges(7,1) * t95 + mrSges(7,2) * t94;
t282 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t152 + t64;
t127 = t175 * t176 - t239;
t120 = t127 * t181;
t66 = t167 * t97 + t172 * t98;
t319 = -t178 * t66 + t120;
t165 = Ifges(6,4) * t257;
t136 = Ifges(7,4) * t146;
t89 = Ifges(7,1) * t147 + Ifges(7,5) * t160 + t136;
t318 = Ifges(6,1) * t258 + Ifges(6,5) * qJD(5) + t180 * t89 + t165;
t281 = qJD(5) * t34;
t8 = t111 * t181 - t178 * t18 - t281;
t141 = qJDD(6) - t151;
t77 = mrSges(7,1) * t141 - mrSges(7,3) * t94;
t78 = -mrSges(7,2) * t141 + mrSges(7,3) * t95;
t316 = -t177 * t77 + t180 * t78;
t314 = -t178 * t8 + t181 * t7;
t313 = m(5) + t334;
t222 = mrSges(6,1) * t181 - mrSges(6,2) * t178;
t312 = -t178 * t245 - t181 * t196 - mrSges(5,1) - t222;
t292 = Ifges(6,4) * t178;
t217 = Ifges(6,2) * t181 + t292;
t311 = Ifges(6,6) * qJD(5) / 0.2e1 + qJD(3) * t217 / 0.2e1 + Ifges(7,5) * t331 + Ifges(7,6) * t332 + Ifges(7,3) * t330;
t209 = t167 * t182 + t172 * t179;
t122 = t209 * t170;
t124 = t209 * t175;
t260 = t182 * t172;
t142 = t167 * t179 - t260;
t70 = t171 * (t124 * t173 - t142 * t168) + t122 * t176;
t220 = mrSges(7,1) * t177 + mrSges(7,2) * t180;
t310 = pkin(9) * t334 - mrSges(5,2) + t220;
t309 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t308 = qJD(3) ^ 2;
t290 = Ifges(7,4) * t147;
t88 = Ifges(7,2) * t146 + Ifges(7,6) * t160 + t290;
t307 = -t88 / 0.2e1;
t306 = t94 / 0.2e1;
t305 = t95 / 0.2e1;
t304 = t141 / 0.2e1;
t302 = t147 / 0.2e1;
t6 = -qJDD(5) * pkin(5) - t8;
t297 = t178 * t6;
t291 = Ifges(6,4) * t181;
t289 = Ifges(7,4) * t177;
t288 = Ifges(7,4) * t180;
t286 = t178 * t56;
t278 = t127 * t178;
t169 = sin(pkin(11));
t174 = cos(pkin(11));
t266 = t174 * t176;
t129 = t168 * t266 + t169 * t173;
t277 = t129 * t179;
t272 = t169 * t176;
t131 = -t168 * t272 + t173 * t174;
t276 = t131 * t179;
t273 = t169 * t171;
t268 = t171 * t174;
t267 = t171 * t175;
t264 = t177 * t178;
t262 = t178 * t180;
t251 = qJD(6) * t178;
t248 = pkin(3) * t265;
t247 = Ifges(7,5) * t94 + Ifges(7,6) * t95 + Ifges(7,3) * t141;
t243 = mrSges(6,3) * t257;
t238 = t171 * t269;
t232 = t177 * t253;
t228 = t249 / 0.2e1;
t227 = m(3) + m(4) + t313;
t225 = mrSges(4,1) * t97 - mrSges(4,2) * t98;
t219 = Ifges(7,1) * t180 - t289;
t218 = Ifges(7,1) * t177 + t288;
t216 = -Ifges(7,2) * t177 + t288;
t215 = Ifges(7,2) * t180 + t289;
t214 = Ifges(6,5) * t181 - Ifges(6,6) * t178;
t213 = Ifges(7,5) * t180 - Ifges(7,6) * t177;
t212 = Ifges(7,5) * t177 + Ifges(7,6) * t180;
t39 = t181 * t66 + t278;
t65 = t167 * t98 - t172 * t97;
t20 = t177 * t65 + t180 * t39;
t19 = -t177 * t39 + t180 * t65;
t33 = t114 * t181 - t178 * t49;
t105 = t122 * t181 + t175 * t178;
t121 = t167 * t270 - t170 * t260;
t76 = t105 * t180 + t121 * t177;
t75 = -t105 * t177 + t121 * t180;
t104 = t122 * t178 - t181 * t175;
t130 = -t168 * t174 - t173 * t272;
t208 = t130 * t248 + (t169 * t238 - t276) * pkin(3);
t128 = -t168 * t169 + t173 * t266;
t206 = -t128 * t175 + t170 * t268;
t205 = t130 * t175 + t169 * t271;
t48 = -qJD(3) * pkin(4) - t52;
t200 = t48 * (mrSges(6,1) * t178 + mrSges(6,2) * t181);
t199 = t97 * pkin(3);
t198 = t178 * (Ifges(6,1) * t181 - t292);
t193 = -t177 * t251 + t180 * t253;
t192 = t178 * t250 + t232;
t40 = t122 * t268 - t124 * t128 + t129 * t142;
t45 = t122 * t273 + t124 * t130 - t131 * t142;
t189 = Ifges(7,5) * t178 + t181 * t219;
t188 = Ifges(7,6) * t178 + t181 * t216;
t187 = Ifges(7,3) * t178 + t181 * t213;
t186 = t128 * t248 + (-t174 * t238 - t277) * pkin(3);
t164 = -pkin(4) - t300;
t156 = -qJD(5) * mrSges(6,2) + t243;
t149 = t224 * qJD(3);
t148 = t222 * qJD(3);
t125 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t151;
t123 = t142 * t175;
t119 = t142 * t170 * qJD(3);
t118 = qJD(3) * t122;
t116 = mrSges(7,1) * t160 - mrSges(7,3) * t147;
t115 = -mrSges(7,2) * t160 + mrSges(7,3) * t146;
t108 = -mrSges(6,1) * t151 + mrSges(6,2) * t152;
t100 = -t130 * t170 + t169 * t267;
t99 = -t128 * t170 - t174 * t267;
t93 = t98 * qJD(3);
t92 = t97 * qJD(3);
t74 = qJD(5) * t105 - t119 * t178;
t73 = -qJD(5) * t104 - t119 * t181;
t69 = -t121 * t176 + (-t123 * t173 - t168 * t209) * t171;
t63 = -t167 * t93 + t172 * t92;
t62 = t167 * t92 + t172 * t93;
t51 = t94 * Ifges(7,1) + t95 * Ifges(7,4) + t141 * Ifges(7,5);
t50 = t94 * Ifges(7,4) + t95 * Ifges(7,2) + t141 * Ifges(7,6);
t47 = t181 * t70 + t278;
t44 = -t121 * t273 - t123 * t130 - t131 * t209;
t41 = t121 * t268 - t123 * t128 - t129 * t209;
t31 = -qJD(5) * pkin(5) - t33;
t30 = t100 * t178 + t181 * t45;
t28 = t178 * t99 - t181 * t40;
t26 = t149 * t177 + t180 * t33;
t25 = t149 * t180 - t177 * t33;
t24 = -qJD(6) * t76 + t118 * t180 - t177 * t73;
t23 = qJD(6) * t75 + t118 * t177 + t180 * t73;
t13 = qJD(5) * t39 + t178 * t63;
t12 = qJD(5) * t319 + t181 * t63;
t4 = -qJD(6) * t20 - t12 * t177 + t180 * t62;
t3 = qJD(6) * t19 + t12 * t180 + t177 * t62;
t14 = [m(2) * qJDD(1) + t65 * t108 + t3 * t115 + t4 * t116 + t12 * t156 + t39 * t125 - t62 * t148 + t19 * t77 + t20 * t78 - t282 * t319 + t259 * t13 + (-mrSges(5,1) * t65 - mrSges(5,2) * t66 + t225) * qJDD(3) + (-mrSges(4,1) * t93 - mrSges(5,1) * t62 - mrSges(4,2) * t92 - mrSges(5,2) * t63) * qJD(3) + (-m(2) - t227) * g(3) + m(7) * (t1 * t20 + t10 * t3 + t13 * t31 + t19 * t2 - t319 * t6 + t4 * t9) + m(6) * (t12 * t34 - t13 * t33 + t17 * t65 + t319 * t8 + t39 * t7 + t48 * t62) + m(5) * (t111 * t127 - t21 * t65 + t22 * t66 - t52 * t62 + t53 * t63) + m(4) * (t117 * t127 + t60 * t98 + t61 * t97 - t85 * t93 + t86 * t92) + m(3) * (t158 * t176 + (t168 ^ 2 + t173 ^ 2) * t171 ^ 2 * qJDD(1)); t105 * t125 + t121 * t108 + t23 * t115 + t24 * t116 - t118 * t148 + t73 * t156 + t75 * t77 + t76 * t78 + t259 * t74 + (-mrSges(4,1) * t179 - mrSges(4,2) * t182) * t308 * t170 + t282 * t104 + (-mrSges(5,1) * t118 + mrSges(5,2) * t119) * qJD(3) + (-t121 * mrSges(5,1) - t122 * mrSges(5,2) + (mrSges(4,1) * t182 - mrSges(4,2) * t179) * t170) * qJDD(3) + m(7) * (t1 * t76 + t10 * t23 + t104 * t6 + t2 * t75 + t24 * t9 + t31 * t74) + m(6) * (-t104 * t8 + t105 * t7 + t118 * t48 + t121 * t17 - t33 * t74 + t34 * t73) + m(3) * t158 + m(5) * (t111 * t175 - t118 * t52 - t119 * t53 - t121 * t21 + t122 * t22) + (t117 * t175 + (t179 * t60 + t182 * t61 + (-t179 * t85 + t182 * t86) * qJD(3)) * t170) * m(4) + (-t176 * g(3) + (-g(1) * t169 + g(2) * t174) * t171) * t227; -(t177 * t89 + t180 * t88) * t251 / 0.2e1 + (-t167 * t280 - t22) * mrSges(5,2) + (t172 * t280 + t21) * mrSges(5,1) + t198 * t228 + t318 * t253 / 0.2e1 + qJD(5) ^ 2 * t214 / 0.2e1 + (t9 * mrSges(7,1) - t10 * mrSges(7,2) - t311) * t255 - t17 * t222 + (t164 * t17 + ((-t178 * t34 - t181 * t33) * qJD(5) + t314) * t163 - t48 * t55 - (-t178 * t33 + t181 * t34) * t56) * m(6) + t282 * t163 * t178 + (t52 * t55 - t53 * t56 + (t167 * t22 + t172 * t21) * pkin(3)) * m(5) + (-g(1) * t45 + g(2) * t40 - g(3) * t70 - t253 * t33 - t255 * t34 + t314) * mrSges(6,3) + (-t1 * t264 - t10 * t192 - t193 * t9 - t2 * t262) * mrSges(7,3) + (Ifges(6,4) * t325 + Ifges(6,2) * t326 + (-Ifges(6,2) * t178 + t291) * t228 - Ifges(7,3) * t304 - Ifges(7,6) * t305 - Ifges(7,5) * t306 - t247 / 0.2e1 - t56 * t156 + t163 * t125 + t309) * t181 + (Ifges(6,1) * t152 + Ifges(6,4) * t326 + t213 * t304 + t216 * t305 + t219 * t306) * t178 - t259 * (-t163 * t253 + t286) + t31 * (mrSges(7,1) * t192 + mrSges(7,2) * t193) + (mrSges(4,1) * t86 + mrSges(5,1) * t55 + mrSges(4,2) * t85 + mrSges(5,2) * t56) * qJD(3) + t323 * t115 + (t1 * t102 + t101 * t2 + (t253 * t31 + t297) * t163 - t286 * t31 + t324 * t9 + t323 * t10) * m(7) + t324 * t116 + (-m(5) * t208 - (t182 * t205 - t276) * mrSges(4,1) - (-t131 * t182 - t179 * t205) * mrSges(4,2) - t334 * (t44 * pkin(4) + t208) - t310 * t45 + t312 * t44) * g(1) + (-m(5) * t199 - t225 - t334 * (t69 * pkin(4) + t199) - t310 * t70 + t312 * t69) * g(3) + (-m(5) * t186 - (-t182 * t206 - t277) * mrSges(4,1) - (-t129 * t182 + t179 * t206) * mrSges(4,2) - t334 * (t41 * pkin(4) + t186) + t310 * t40 + t312 * t41) * g(2) + qJDD(5) * (Ifges(6,5) * t178 + Ifges(6,6) * t181) + t164 * t108 + t55 * t148 - t50 * t264 / 0.2e1 + t51 * t262 / 0.2e1 + qJD(5) * t200 + t220 * t297 + t291 * t325 + t217 * t326 + (qJD(5) * t189 - t218 * t251) * t302 + t232 * t307 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) - t156 * t234 - t60 * mrSges(4,2) + t61 * mrSges(4,1) + t146 * (qJD(5) * t188 - t215 * t251) / 0.2e1 + t160 * (qJD(5) * t187 - t212 * t251) / 0.2e1 + t101 * t77 + t102 * t78; t111 * m(5) + ((t115 * t180 - t116 * t177 + t156) * qJD(5) + m(6) * (t8 + t281) + m(7) * (t10 * t254 - t256 * t9 - t6) - t282) * t181 + (t125 + (-t177 * t115 - t180 * t116) * qJD(6) + t259 * qJD(5) + m(6) * (-qJD(5) * t33 + t7) + m(7) * (qJD(5) * t31 + t329) + t316) * t178 + t313 * (-g(1) * t100 - g(2) * t99 - g(3) * t127); (t146 * t216 + t147 * t219 + t160 * t213) * qJD(6) / 0.2e1 - (t146 * t188 + t147 * t189 + t160 * t187) * qJD(3) / 0.2e1 + (t257 * t88 + t51) * t177 / 0.2e1 + (t321 * t30 + t322 * (t100 * t181 - t178 * t45)) * g(1) + (t321 * t47 + t322 * (-t178 * t70 + t120)) * g(3) + (-pkin(5) * t6 - t10 * t26 - t25 * t9) * m(7) - (-Ifges(6,2) * t258 + t165 + t318) * t257 / 0.2e1 + (-m(7) * t31 + t244 - t259) * t34 + t311 * t258 + (t321 * t28 + t322 * (t178 * t40 + t181 * t99)) * g(2) + t6 * t221 + (m(7) * ((-t10 * t177 - t180 * t9) * qJD(6) + t223) - t116 * t250 - t115 * t252 + t316) * pkin(10) + t160 * t31 * t220 + (-t200 - t9 * (mrSges(7,1) * t178 - mrSges(7,3) * t261) - t10 * (-mrSges(7,2) * t178 - mrSges(7,3) * t263)) * qJD(3) + t329 * mrSges(7,3) + Ifges(6,3) * qJDD(5) + (t243 - t156) * t33 - t308 * t198 / 0.2e1 + t180 * t50 / 0.2e1 + Ifges(6,5) * t152 + Ifges(6,6) * t151 - t25 * t116 + t212 * t304 + t215 * t305 + t218 * t306 + t252 * t307 - t7 * mrSges(6,2) + t8 * mrSges(6,1) - pkin(5) * t64 - t214 * t249 / 0.2e1 + t89 * t250 / 0.2e1 - t26 * t115; -t31 * (mrSges(7,1) * t147 + mrSges(7,2) * t146) + (Ifges(7,1) * t146 - t290) * t331 + t88 * t302 + (Ifges(7,5) * t146 - Ifges(7,6) * t147) * t330 - t9 * t115 + t10 * t116 - g(1) * ((-t177 * t30 - t180 * t44) * mrSges(7,1) + (t177 * t44 - t180 * t30) * mrSges(7,2)) - g(2) * ((-t177 * t28 - t180 * t41) * mrSges(7,1) + (t177 * t41 - t180 * t28) * mrSges(7,2)) - g(3) * ((-t177 * t47 - t180 * t69) * mrSges(7,1) + (t177 * t69 - t180 * t47) * mrSges(7,2)) + (t10 * t147 + t146 * t9) * mrSges(7,3) + t247 + (-Ifges(7,2) * t147 + t136 + t89) * t332 - t309;];
tau  = t14;
