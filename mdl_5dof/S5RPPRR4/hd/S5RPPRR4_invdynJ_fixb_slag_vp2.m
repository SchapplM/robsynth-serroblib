% Calculate vector of inverse dynamics joint torques for
% S5RPPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% m [6x1]
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
% Datum: 2022-01-23 09:17
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:16:17
% EndTime: 2022-01-23 09:16:34
% DurationCPUTime: 10.12s
% Computational Cost: add. (4551->426), mult. (11546->585), div. (0->0), fcn. (8486->14), ass. (0->209)
t182 = sin(pkin(8));
t184 = cos(pkin(8));
t229 = t182 ^ 2 + t184 ^ 2;
t220 = qJD(1) * qJD(2);
t159 = qJDD(1) * qJ(2) + t220;
t180 = pkin(9) + qJ(4);
t169 = sin(t180);
t255 = pkin(4) * t169;
t181 = sin(pkin(9));
t257 = pkin(3) * t181;
t183 = cos(pkin(9));
t278 = mrSges(4,1) * t181 + mrSges(4,2) * t183;
t291 = -m(6) - m(3);
t294 = (-m(4) + t291) * qJ(2) - m(5) * (qJ(2) + t257) - m(6) * (t255 + t257) + mrSges(2,2) - mrSges(3,3) - t278;
t186 = sin(qJ(4));
t189 = cos(qJ(4));
t141 = t181 * t189 + t183 * t186;
t133 = t141 * qJD(4);
t194 = qJD(1) * t141;
t271 = t184 * t194 - t133;
t234 = t183 * t189;
t199 = t181 * t186 - t234;
t227 = qJD(1) * t184;
t269 = t199 * qJD(4);
t270 = t199 * t227 - t269;
t146 = pkin(2) * t184 + qJ(3) * t182 + pkin(1);
t202 = mrSges(3,1) * t184 - mrSges(3,2) * t182;
t252 = pkin(3) * t183 + pkin(2);
t253 = qJ(3) + pkin(6);
t292 = m(4) * t146 + m(5) * (t184 * t252 + pkin(1)) + mrSges(2,1) + t202 + (m(5) * t253 + mrSges(4,3) + mrSges(5,3) + mrSges(6,3)) * t182;
t68 = (-qJD(1) * t133 - qJDD(1) * t199) * t182;
t290 = Ifges(5,5) * t68;
t185 = sin(qJ(5));
t188 = cos(qJ(5));
t107 = t182 * t194;
t228 = qJD(1) * t182;
t210 = t181 * t228;
t109 = -t186 * t210 + t228 * t234;
t204 = -t107 * t188 - t109 * t185;
t69 = (qJD(1) * t269 - qJDD(1) * t141) * t182;
t21 = qJD(5) * t204 + t185 * t69 + t188 * t68;
t289 = Ifges(6,5) * t21;
t288 = Ifges(5,6) * t69;
t60 = -t107 * t185 + t109 * t188;
t22 = -qJD(5) * t60 - t185 * t68 + t188 * t69;
t287 = Ifges(6,6) * t22;
t217 = qJDD(1) * t184;
t160 = qJDD(4) - t217;
t285 = Ifges(5,3) * t160;
t153 = qJDD(5) + t160;
t284 = Ifges(6,3) * t153;
t195 = t278 * t182;
t283 = -mrSges(5,1) * t107 - mrSges(5,2) * t109 - qJD(1) * t195;
t126 = -qJD(1) * t146 + qJD(2);
t114 = t183 * t126;
t238 = t182 * t183;
t216 = pkin(6) * t238;
t192 = -t216 + (-qJ(2) * t181 - pkin(3)) * t184;
t70 = qJD(1) * t192 + t114;
t211 = qJ(2) * t227;
t88 = t126 * t181 + t183 * t211;
t77 = -pkin(6) * t210 + t88;
t38 = t186 * t70 + t189 * t77;
t239 = t181 * t184;
t224 = qJD(3) * t182;
t94 = -qJD(1) * t224 - qJDD(1) * t146 + qJDD(2);
t73 = -t159 * t239 + t183 * t94;
t48 = (-pkin(3) * t184 - t216) * qJDD(1) + t73;
t218 = qJDD(1) * t182;
t206 = t181 * t218;
t235 = t183 * t184;
t74 = t159 * t235 + t181 * t94;
t55 = -pkin(6) * t206 + t74;
t13 = -qJD(4) * t38 - t186 * t55 + t189 * t48;
t6 = pkin(4) * t160 - pkin(7) * t68 + t13;
t222 = qJD(4) * t189;
t223 = qJD(4) * t186;
t12 = t186 * t48 + t189 * t55 + t222 * t70 - t223 * t77;
t7 = pkin(7) * t69 + t12;
t29 = -pkin(7) * t107 + t38;
t243 = t185 * t29;
t162 = qJD(4) - t227;
t37 = -t186 * t77 + t189 * t70;
t28 = -pkin(7) * t109 + t37;
t25 = pkin(4) * t162 + t28;
t8 = t188 * t25 - t243;
t2 = qJD(5) * t8 + t185 * t6 + t188 * t7;
t242 = t188 * t29;
t9 = t185 * t25 + t242;
t3 = -qJD(5) * t9 - t185 * t7 + t188 * t6;
t281 = -t3 * mrSges(6,1) + t2 * mrSges(6,2);
t280 = -t13 * mrSges(5,1) + t12 * mrSges(5,2);
t279 = mrSges(4,1) * t183 - mrSges(4,2) * t181;
t83 = -t141 * t185 - t188 * t199;
t275 = qJD(5) * t83 + t185 * t271 + t188 * t270;
t84 = t141 * t188 - t185 * t199;
t274 = -qJD(5) * t84 - t185 * t270 + t188 * t271;
t268 = m(6) * pkin(4);
t273 = mrSges(5,1) + t268;
t272 = qJD(1) ^ 2 * t229;
t267 = -t204 / 0.2e1;
t266 = -t60 / 0.2e1;
t265 = t60 / 0.2e1;
t264 = t8 * mrSges(6,3);
t263 = t9 * mrSges(6,3);
t261 = t109 / 0.2e1;
t156 = qJD(5) + t162;
t260 = -t156 / 0.2e1;
t258 = Ifges(6,4) * t60;
t256 = pkin(4) * t109;
t254 = g(3) * t182;
t137 = t183 * t146;
t80 = -t137 + t192;
t101 = qJ(2) * t235 - t146 * t181;
t240 = t181 * t182;
t89 = -pkin(6) * t240 + t101;
t42 = t186 * t80 + t189 * t89;
t172 = qJ(5) + t180;
t166 = sin(t172);
t167 = cos(t172);
t190 = cos(qJ(1));
t187 = sin(qJ(1));
t233 = t184 * t187;
t103 = t166 * t233 + t167 * t190;
t104 = t166 * t190 - t167 * t233;
t251 = -mrSges(6,1) * t103 + mrSges(6,2) * t104;
t232 = t184 * t190;
t105 = -t166 * t232 + t167 * t187;
t106 = t166 * t187 + t167 * t232;
t250 = mrSges(6,1) * t105 - mrSges(6,2) * t106;
t249 = mrSges(5,3) * t107;
t248 = mrSges(5,3) * t109;
t247 = Ifges(5,4) * t109;
t171 = t182 * qJ(2);
t205 = t183 * t218;
t231 = mrSges(4,1) * t206 + mrSges(4,2) * t205;
t142 = pkin(3) * t240 + t171;
t226 = qJD(2) * t182;
t225 = qJD(2) * t184;
t157 = qJ(2) * t228 + qJD(3);
t221 = m(4) + m(5) + m(6);
t139 = t159 * t182 + qJDD(3);
t215 = t284 + t287 + t289;
t214 = t285 + t288 + t290;
t127 = pkin(3) * t210 + t157;
t208 = -t69 * mrSges(5,1) + mrSges(5,2) * t68;
t207 = -t22 * mrSges(6,1) + mrSges(6,2) * t21;
t102 = pkin(3) * t206 + t139;
t41 = -t186 * t89 + t189 * t80;
t203 = -mrSges(3,1) * t217 + mrSges(3,2) * t218;
t201 = -mrSges(6,1) * t166 - mrSges(6,2) * t167;
t123 = t199 * t182;
t35 = -pkin(4) * t184 + pkin(7) * t123 + t41;
t122 = t141 * t182;
t36 = -pkin(7) * t122 + t42;
t14 = -t185 * t36 + t188 * t35;
t15 = t185 * t35 + t188 * t36;
t71 = -t122 * t188 + t123 * t185;
t72 = -t122 * t185 - t123 * t188;
t170 = cos(t180);
t200 = (pkin(4) * t170 + t252) * t184 - (-pkin(7) - t253) * t182;
t198 = t215 - t281;
t197 = -mrSges(4,1) * t184 - mrSges(4,3) * t238;
t196 = mrSges(4,2) * t184 - mrSges(4,3) * t240;
t117 = -t169 * t232 + t170 * t187;
t115 = t169 * t233 + t170 * t190;
t128 = -t181 * t225 - t183 * t224;
t129 = -t181 * t224 + t183 * t225;
t30 = t128 * t186 + t129 * t189 + t222 * t80 - t223 * t89;
t31 = -qJD(4) * t42 + t128 * t189 - t129 * t186;
t168 = -qJDD(1) * pkin(1) + qJDD(2);
t135 = t197 * qJD(1);
t134 = t196 * qJD(1);
t131 = t197 * qJDD(1);
t130 = t196 * qJDD(1);
t118 = t169 * t187 + t170 * t232;
t116 = t169 * t190 - t170 * t233;
t112 = t182 * t269;
t111 = t182 * t133;
t100 = -qJ(2) * t239 - t137;
t99 = Ifges(5,4) * t107;
t90 = -pkin(4) * t112 + t226;
t87 = -t181 * t211 + t114;
t86 = mrSges(5,1) * t162 - t248;
t85 = -mrSges(5,2) * t162 - t249;
t82 = pkin(4) * t122 + t142;
t75 = pkin(4) * t107 + t127;
t54 = Ifges(6,4) * t204;
t52 = Ifges(5,1) * t109 + Ifges(5,5) * t162 - t99;
t51 = -Ifges(5,2) * t107 + Ifges(5,6) * t162 + t247;
t50 = -mrSges(5,2) * t160 + mrSges(5,3) * t69;
t49 = mrSges(5,1) * t160 - mrSges(5,3) * t68;
t45 = mrSges(6,1) * t156 - mrSges(6,3) * t60;
t44 = -mrSges(6,2) * t156 + mrSges(6,3) * t204;
t43 = -pkin(4) * t69 + t102;
t34 = -qJD(5) * t72 + t111 * t185 + t112 * t188;
t33 = qJD(5) * t71 - t111 * t188 + t112 * t185;
t32 = -mrSges(6,1) * t204 + mrSges(6,2) * t60;
t27 = Ifges(6,1) * t60 + Ifges(6,5) * t156 + t54;
t26 = Ifges(6,2) * t204 + Ifges(6,6) * t156 + t258;
t24 = pkin(7) * t111 + t31;
t23 = pkin(7) * t112 + t30;
t17 = -mrSges(6,2) * t153 + mrSges(6,3) * t22;
t16 = mrSges(6,1) * t153 - mrSges(6,3) * t21;
t11 = t188 * t28 - t243;
t10 = -t185 * t28 - t242;
t5 = -qJD(5) * t15 - t185 * t23 + t188 * t24;
t4 = qJD(5) * t14 + t185 * t24 + t188 * t23;
t1 = [(-mrSges(5,2) * t102 + mrSges(5,3) * t13 - Ifges(5,1) * t68 - Ifges(5,4) * t69 - Ifges(5,5) * t160) * t123 + t231 * t171 + t139 * t195 + t73 * t197 + t74 * t196 + (mrSges(6,2) * t43 - mrSges(6,3) * t3 + Ifges(6,1) * t21 + Ifges(6,4) * t22 + Ifges(6,5) * t153) * t72 + t34 * t263 + (Ifges(6,1) * t33 + Ifges(6,4) * t34) * t265 + (Ifges(3,4) * t218 - Ifges(4,5) * t205 + Ifges(4,6) * t206 - t289 / 0.2e1 - t287 / 0.2e1 - t290 / 0.2e1 - t288 / 0.2e1 - t284 / 0.2e1 - t285 / 0.2e1 + (Ifges(3,2) + Ifges(4,3)) * t217 + t280 + t281) * t184 + m(5) * (t102 * t142 + t12 * t42 + t13 * t41 + t30 * t38 + t31 * t37) + (m(5) * t127 - t283) * t226 + (-t116 * mrSges(5,1) - t104 * mrSges(6,1) - t115 * mrSges(5,2) - t103 * mrSges(6,2) + t279 * t233 + t294 * t190 + (-m(6) * (-pkin(1) - t200) + m(3) * pkin(1) + t292) * t187) * g(1) + (-t118 * mrSges(5,1) - t106 * mrSges(6,1) - t117 * mrSges(5,2) - t105 * mrSges(6,2) - t279 * t232 + (-m(6) * t200 + pkin(1) * t291 - t292) * t190 + t294 * t187) * g(2) + m(3) * (-pkin(1) * t168 + (t159 + t220) * qJ(2) * t229) + t204 * (Ifges(6,4) * t33 + Ifges(6,2) * t34) / 0.2e1 - (t215 + t214) * t184 / 0.2e1 + (-Ifges(5,1) * t111 + Ifges(5,4) * t112) * t261 + (t111 * t37 + t112 * t38) * mrSges(5,3) - t107 * (-Ifges(5,4) * t111 + Ifges(5,2) * t112) / 0.2e1 + t127 * (-mrSges(5,1) * t112 - mrSges(5,2) * t111) + t162 * (-Ifges(5,5) * t111 + Ifges(5,6) * t112) / 0.2e1 + (-mrSges(6,1) * t43 + mrSges(6,3) * t2 + Ifges(6,4) * t21 + Ifges(6,2) * t22 + Ifges(6,6) * t153) * t71 + m(6) * (t14 * t3 + t15 * t2 + t4 * t9 + t43 * t82 + t5 * t8 + t75 * t90) + (mrSges(5,1) * t102 - mrSges(5,3) * t12 - Ifges(5,4) * t68 - Ifges(5,2) * t69 - Ifges(5,6) * t160) * t122 - t168 * t202 - pkin(1) * t203 - t33 * t264 + t82 * t207 + t142 * t208 + m(4) * (t100 * t73 + t101 * t74 + t128 * t87 + t129 * t88) + 0.2e1 * t229 * t159 * mrSges(3,3) + Ifges(2,3) * qJDD(1) + t14 * t16 + t15 * t17 + t33 * t27 / 0.2e1 + t34 * t26 / 0.2e1 + (m(4) * (qJ(2) * t139 + qJD(2) * t157) + Ifges(3,1) * t218 + (Ifges(4,1) * t183 - Ifges(4,4) * t181) * t205 - (Ifges(4,4) * t183 - Ifges(4,2) * t181) * t206 + (-Ifges(4,5) * t183 + Ifges(4,6) * t181 + Ifges(3,4)) * t217) * t182 + t4 * t44 + t5 * t45 + t41 * t49 + t42 * t50 + t75 * (-mrSges(6,1) * t34 + mrSges(6,2) * t33) + t30 * t85 + t31 * t86 + t90 * t32 - t111 * t52 / 0.2e1 + t112 * t51 / 0.2e1 + t101 * t130 + t100 * t131 + t129 * t134 + t128 * t135 + t156 * (Ifges(6,5) * t33 + Ifges(6,6) * t34) / 0.2e1; t271 * t86 + t270 * t85 + t274 * t45 + t275 * t44 + ((-t134 * t183 + t135 * t181) * t184 + (-t32 + t283) * t182) * qJD(1) + t203 - mrSges(3,3) * t272 + t83 * t16 + t84 * t17 - t199 * t49 + t141 * t50 + t181 * t130 + t183 * t131 + (-g(1) * t187 + g(2) * t190) * (m(3) + t221) + (t2 * t84 - t228 * t75 + t274 * t8 + t275 * t9 + t3 * t83) * m(6) + (t12 * t141 - t127 * t228 - t13 * t199 + t270 * t38 + t271 * t37) * m(5) + (t74 * t181 + t73 * t183 - (t157 * t182 + (-t181 * t87 + t183 * t88) * t184) * qJD(1)) * m(4) + (-qJ(2) * t272 + t168) * m(3); t107 * t85 + t109 * t86 - t204 * t44 + t60 * t45 + t221 * t184 * g(3) + m(4) * t139 + ((t183 * t135 + t181 * t134 - m(4) * (-t181 * t88 - t183 * t87)) * qJD(1) + (-g(1) * t190 - g(2) * t187) * t221) * t182 + t207 + t208 + t231 + (-t204 * t9 + t60 * t8 + t43) * m(6) + (t107 * t38 + t109 * t37 + t102) * m(5); (-mrSges(5,2) * t116 + t115 * t273 - t251) * g(2) + (mrSges(5,2) * t118 - t117 * t273 - t250) * g(1) + t214 - t32 * t256 - m(6) * (t10 * t8 + t11 * t9 + t256 * t75) + (m(6) * t255 + mrSges(5,1) * t169 + mrSges(5,2) * t170 - t201) * t254 + t51 * t261 + (-Ifges(5,2) * t109 + t52 - t99) * t107 / 0.2e1 + (Ifges(6,5) * t260 + t264 + Ifges(6,1) * t266 + Ifges(6,4) * t267 - t27 / 0.2e1 - t75 * mrSges(6,2)) * t204 - (Ifges(6,6) * t260 - t263 + Ifges(6,4) * t266 + Ifges(6,2) * t267 - t26 / 0.2e1 + t75 * mrSges(6,1)) * t60 - t109 * (-Ifges(5,1) * t107 - t247) / 0.2e1 - t127 * (mrSges(5,1) * t109 - mrSges(5,2) * t107) - t162 * (-Ifges(5,5) * t107 - Ifges(5,6) * t109) / 0.2e1 + (-t249 - t85) * t37 + t198 - t280 + ((-t185 * t45 + t188 * t44) * qJD(5) + t188 * t16 + t185 * t17) * pkin(4) + (t185 * t2 + t188 * t3 + (-t185 * t8 + t188 * t9) * qJD(5)) * t268 - t11 * t44 - t10 * t45 + (t248 + t86) * t38; -t75 * (mrSges(6,1) * t60 + mrSges(6,2) * t204) + (Ifges(6,1) * t204 - t258) * t266 + t26 * t265 + (Ifges(6,5) * t204 - Ifges(6,6) * t60) * t260 - t8 * t44 + t9 * t45 - g(1) * t250 - g(2) * t251 - t201 * t254 + (t204 * t8 + t60 * t9) * mrSges(6,3) + t198 + (-Ifges(6,2) * t60 + t27 + t54) * t267;];
tau = t1;
