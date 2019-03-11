% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR5
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
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
% Datum: 2019-03-09 01:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR5_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:37:13
% EndTime: 2019-03-09 01:37:22
% DurationCPUTime: 5.80s
% Computational Cost: add. (3362->441), mult. (5423->598), div. (0->0), fcn. (2765->8), ass. (0->202)
t124 = sin(qJ(6));
t127 = cos(qJ(6));
t197 = qJD(5) * t127;
t125 = sin(qJ(5));
t201 = qJD(1) * t125;
t81 = -t124 * t201 + t197;
t279 = -t81 / 0.2e1;
t199 = qJD(5) * t124;
t82 = t127 * t201 + t199;
t278 = -t82 / 0.2e1;
t128 = cos(qJ(5));
t200 = qJD(1) * t128;
t101 = qJD(6) - t200;
t277 = -t101 / 0.2e1;
t120 = sin(pkin(9));
t174 = t125 * t197;
t121 = cos(pkin(9));
t205 = t124 * t128;
t63 = -t120 * t205 - t121 * t127;
t203 = t127 * t128;
t66 = t120 * t124 + t121 * t203;
t276 = t66 * qJD(1) + qJD(6) * t63 - t120 * t174;
t146 = -t120 * t127 + t121 * t205;
t65 = t120 * t203 - t121 * t124;
t275 = -t65 * qJD(1) - qJD(6) * t146 - t121 * t174;
t198 = qJD(5) * t125;
t176 = t124 * t198;
t274 = -t146 * qJD(1) - qJD(6) * t65 + t120 * t176;
t273 = -t63 * qJD(1) - qJD(6) * t66 + t121 * t176;
t182 = mrSges(6,3) * t201;
t225 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t81 + mrSges(7,2) * t82 + t182;
t272 = t125 * t225;
t188 = qJD(1) * qJD(5);
t91 = qJDD(1) * t125 + t128 * t188;
t36 = qJD(6) * t81 + qJDD(5) * t124 + t127 * t91;
t37 = -qJD(6) * t82 + qJDD(5) * t127 - t124 * t91;
t13 = -mrSges(7,1) * t37 + mrSges(7,2) * t36;
t226 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t91 + t13;
t271 = t226 * t125;
t108 = qJD(1) * qJ(2) + qJD(3);
t98 = qJD(1) * pkin(3) + t108;
t123 = -pkin(1) - qJ(3);
t99 = qJD(1) * t123 + qJD(2);
t49 = t120 * t98 + t121 * t99;
t46 = qJD(1) * pkin(7) + t49;
t33 = qJD(4) * t125 + t128 * t46;
t22 = qJD(5) * pkin(8) + t33;
t166 = pkin(5) * t128 + pkin(8) * t125;
t150 = -pkin(4) - t166;
t48 = -t120 * t99 + t121 * t98;
t24 = qJD(1) * t150 - t48;
t10 = t124 * t24 + t127 * t22;
t118 = qJD(1) * qJD(2);
t257 = qJDD(1) * qJ(2) + t118;
t97 = qJDD(3) + t257;
t83 = qJDD(1) * pkin(3) + t97;
t189 = qJD(1) * qJD(3);
t84 = qJDD(1) * t123 + qJDD(2) - t189;
t38 = -t120 * t84 + t121 * t83;
t34 = -qJDD(1) * pkin(4) - t38;
t90 = qJDD(1) * t128 - t125 * t188;
t14 = -pkin(5) * t90 - pkin(8) * t91 + t34;
t196 = qJD(5) * t128;
t39 = t120 * t83 + t121 * t84;
t35 = qJDD(1) * pkin(7) + t39;
t11 = qJD(4) * t196 + t125 * qJDD(4) + t128 * t35 - t198 * t46;
t7 = qJDD(5) * pkin(8) + t11;
t9 = -t124 * t22 + t127 * t24;
t1 = qJD(6) * t9 + t124 * t14 + t127 * t7;
t2 = -qJD(6) * t10 - t124 * t7 + t127 * t14;
t164 = t1 * t127 - t124 * t2;
t193 = qJD(6) * t127;
t195 = qJD(6) * t124;
t270 = -t10 * t195 - t9 * t193 + t164;
t241 = -m(6) / 0.2e1;
t242 = -m(5) / 0.2e1;
t32 = qJD(4) * t128 - t125 * t46;
t255 = -t125 * t32 + t128 * t33;
t21 = -qJD(5) * pkin(5) - t32;
t263 = m(7) * t21;
t268 = t241 * t255 - t125 * t263 / 0.2e1 + t49 * t242;
t130 = qJD(1) ^ 2;
t192 = qJDD(1) * mrSges(5,1);
t181 = mrSges(6,3) * t200;
t95 = -qJD(5) * mrSges(6,2) + t181;
t209 = t128 * t95;
t47 = -mrSges(6,1) * t90 + mrSges(6,2) * t91;
t267 = -t130 * mrSges(5,2) + (t209 + t272) * qJD(1) + t192 - t47;
t266 = t90 / 0.2e1;
t265 = t91 / 0.2e1;
t264 = -m(7) - m(6);
t262 = -mrSges(4,1) - mrSges(3,3);
t261 = mrSges(3,2) - mrSges(4,3);
t260 = mrSges(6,3) - mrSges(5,2);
t163 = mrSges(6,1) * t128 - mrSges(6,2) * t125;
t259 = -mrSges(5,1) - t163;
t107 = Ifges(6,4) * t200;
t75 = Ifges(7,4) * t81;
t27 = Ifges(7,1) * t82 + Ifges(7,5) * t101 + t75;
t258 = Ifges(6,1) * t201 + Ifges(6,5) * qJD(5) + t127 * t27 + t107;
t208 = qJD(5) * t33;
t12 = qJDD(4) * t128 - t125 * t35 - t208;
t78 = qJDD(6) - t90;
t17 = mrSges(7,1) * t78 - mrSges(7,3) * t36;
t18 = -mrSges(7,2) * t78 + mrSges(7,3) * t37;
t254 = -t124 * t17 + t127 * t18;
t253 = t11 * t128 - t12 * t125;
t252 = -mrSges(2,1) + t261;
t251 = mrSges(2,2) + t262;
t223 = Ifges(6,4) * t125;
t157 = t128 * Ifges(6,2) + t223;
t250 = Ifges(6,6) * qJD(5) / 0.2e1 + qJD(1) * t157 / 0.2e1 + Ifges(7,5) * t278 + Ifges(7,6) * t279 + Ifges(7,3) * t277;
t248 = t225 + t263;
t246 = -t2 * mrSges(7,1) + t1 * mrSges(7,2);
t45 = -qJD(1) * pkin(4) - t48;
t245 = -t241 * t45 + t242 * t48;
t243 = 0.2e1 * qJD(1);
t231 = Ifges(7,4) * t82;
t26 = Ifges(7,2) * t81 + Ifges(7,6) * t101 + t231;
t240 = -t26 / 0.2e1;
t239 = t36 / 0.2e1;
t238 = t37 / 0.2e1;
t237 = t78 / 0.2e1;
t235 = t82 / 0.2e1;
t224 = m(4) * qJD(1);
t222 = Ifges(6,4) * t128;
t221 = Ifges(7,4) * t124;
t220 = Ifges(7,4) * t127;
t122 = qJ(2) + pkin(3);
t77 = t120 * t122 + t121 * t123;
t74 = pkin(7) + t77;
t210 = t128 * t74;
t207 = qJDD(1) * pkin(1);
t206 = t124 * t125;
t204 = t125 * t127;
t126 = sin(qJ(1));
t129 = cos(qJ(1));
t202 = t129 * pkin(1) + t126 * qJ(2);
t194 = qJD(6) * t125;
t191 = qJDD(1) * mrSges(5,2);
t190 = -m(5) + t264;
t186 = Ifges(7,5) * t36 + Ifges(7,6) * t37 + Ifges(7,3) * t78;
t184 = -m(4) + t190;
t183 = m(7) * pkin(8) + mrSges(7,3);
t180 = t74 * t198;
t177 = t129 * qJ(3) + t202;
t175 = t124 * t196;
t170 = t188 / 0.2e1;
t76 = -t120 * t123 + t121 * t122;
t168 = t126 * pkin(3) + t177;
t165 = pkin(5) * t125 - pkin(8) * t128;
t92 = qJD(2) * t121 + qJD(3) * t120;
t167 = qJD(5) * t165 - qJD(6) * t210 - t92;
t162 = mrSges(6,1) * t125 + mrSges(6,2) * t128;
t161 = -mrSges(7,1) * t127 + mrSges(7,2) * t124;
t160 = t124 * mrSges(7,1) + t127 * mrSges(7,2);
t159 = Ifges(7,1) * t127 - t221;
t158 = Ifges(7,1) * t124 + t220;
t156 = -Ifges(7,2) * t124 + t220;
t155 = Ifges(7,2) * t127 + t221;
t154 = Ifges(6,5) * t128 - Ifges(6,6) * t125;
t153 = Ifges(7,5) * t127 - Ifges(7,6) * t124;
t152 = Ifges(7,5) * t124 + Ifges(7,6) * t127;
t113 = t129 * qJ(2);
t149 = t129 * pkin(3) + t123 * t126 + t113;
t8 = -qJDD(5) * pkin(5) - t12;
t148 = t125 * t8 + t196 * t21;
t144 = t45 * t162;
t143 = t125 * (Ifges(6,1) * t128 - t223);
t141 = m(7) * pkin(5) - t161;
t139 = -t124 * t194 + t127 * t196;
t138 = t125 * t193 + t175;
t52 = t150 - t76;
t93 = qJD(2) * t120 - qJD(3) * t121;
t137 = -qJD(6) * t52 - t128 * t93 + t180;
t136 = Ifges(7,5) * t125 + t128 * t159;
t135 = Ifges(7,6) * t125 + t128 * t156;
t134 = Ifges(7,3) * t125 + t128 * t153;
t133 = (-t125 * t33 - t128 * t32) * qJD(5) + t253;
t132 = t125 * t183 + t128 * t141;
t58 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t90;
t85 = t163 * qJD(1);
t131 = -t130 * mrSges(5,1) - t191 - qJD(1) * t85 + t128 * t58 + t271 + (-t125 * t95 + t225 * t128) * qJD(5);
t109 = qJDD(2) - t207;
t89 = t165 * qJD(1);
t80 = t120 * t129 + t121 * t126;
t79 = -t126 * t120 + t121 * t129;
t73 = -pkin(4) - t76;
t72 = t160 * t125;
t51 = mrSges(7,1) * t101 - mrSges(7,3) * t82;
t50 = -mrSges(7,2) * t101 + mrSges(7,3) * t81;
t31 = -t124 * t79 + t203 * t80;
t30 = -t127 * t79 - t205 * t80;
t20 = t124 * t52 + t203 * t74;
t19 = t127 * t52 - t205 * t74;
t16 = t124 * t89 + t127 * t32;
t15 = -t124 * t32 + t127 * t89;
t6 = t36 * Ifges(7,1) + t37 * Ifges(7,4) + t78 * Ifges(7,5);
t5 = t36 * Ifges(7,4) + t37 * Ifges(7,2) + t78 * Ifges(7,6);
t4 = t124 * t137 + t127 * t167;
t3 = t124 * t167 - t127 * t137;
t23 = [(-t1 * t206 - t10 * t138 - t139 * t9 - t2 * t204) * mrSges(7,3) + (-qJD(1) * t93 - t39) * mrSges(5,2) + (t109 - t207) * mrSges(3,2) - (t124 * t27 + t127 * t26) * t194 / 0.2e1 + t143 * t170 + qJD(5) ^ 2 * t154 / 0.2e1 + (qJD(1) * t92 + t38) * mrSges(5,1) + qJDD(5) * (Ifges(6,5) * t125 + Ifges(6,6) * t128) + t92 * t85 + 0.2e1 * t257 * mrSges(3,3) + m(4) * (qJ(2) * t97 + qJD(2) * t108 - qJD(3) * t99 + t123 * t84) + m(5) * (t38 * t76 + t39 * t77 + t48 * t92 + t49 * t93) + qJD(5) * t144 + t222 * t265 + t157 * t266 + t74 * t271 + t93 * t272 + m(7) * (t1 * t20 + t10 * t3 + t19 * t2 + t4 * t9) + (-mrSges(4,3) * t123 + Ifges(3,1) + Ifges(4,2) + Ifges(2,3) + Ifges(5,3)) * qJDD(1) + t8 * t72 + t73 * t47 + t3 * t50 + t4 * t51 + t19 * t17 + t20 * t18 + t6 * t204 / 0.2e1 + (-t84 + t189) * mrSges(4,3) + (-t198 * t33 + t253) * mrSges(6,3) + (t258 / 0.2e1 - t32 * mrSges(6,3) + t248 * t74) * t196 - t34 * t163 + t21 * (mrSges(7,1) * t138 + mrSges(7,2) * t139) + t81 * (qJD(5) * t135 - t155 * t194) / 0.2e1 + t101 * (qJD(5) * t134 - t152 * t194) / 0.2e1 - t77 * t191 + m(3) * (-pkin(1) * t109 + (t257 + t118) * qJ(2)) - t95 * t180 - t5 * t206 / 0.2e1 + (-m(5) * t149 + t264 * (t79 * pkin(4) + t149) + (-m(3) - m(4)) * t113 + (pkin(7) * t264 - t160 - t260) * t80 + (-t132 + t259) * t79 + t251 * t129 + (m(3) * pkin(1) - m(4) * t123 - t252) * t126) * g(1) + (-m(3) * t202 - m(4) * t177 - m(5) * t168 - t31 * mrSges(7,1) - t30 * mrSges(7,2) + t260 * t79 + t264 * (t80 * pkin(4) - pkin(7) * t79 + t168) + (-m(7) * t166 - t125 * mrSges(7,3) + t259) * t80 + t252 * t129 + t251 * t126) * g(2) + (Ifges(6,4) * t265 + Ifges(6,2) * t266 + (-Ifges(6,2) * t125 + t222) * t170 - t186 / 0.2e1 - Ifges(7,3) * t237 - Ifges(7,6) * t238 - Ifges(7,5) * t239 + t246) * t128 + (Ifges(6,1) * t91 + Ifges(6,4) * t266 + t153 * t237 + t156 * t238 + t159 * t239 + m(7) * (t21 * t93 + t74 * t8)) * t125 + (t97 + t257) * mrSges(4,1) + m(6) * (t133 * t74 + t255 * t93 + t34 * t73 - t45 * t92) + (t9 * mrSges(7,1) - t10 * mrSges(7,2) - t250) * t198 + t76 * t192 + t93 * t209 + t58 * t210 + (qJD(5) * t136 - t158 * t194) * t235 + t175 * t240; -t108 * t224 - t146 * t17 + t66 * t18 + t273 * t51 + t275 * t50 + t261 * qJDD(1) + (-m(3) * qJ(2) + t262) * t130 - t267 * t120 + t131 * t121 + m(6) * (t120 * t34 + t121 * t133) + m(5) * (-t120 * t38 + t121 * t39) + m(3) * t109 + m(4) * t84 + (t120 * t268 + t245 * t121) * t243 + (-t126 * g(1) + t129 * g(2)) * (m(3) - t184) + (t1 * t66 + t10 * t275 + t121 * t148 - t2 * t146 + t273 * t9) * m(7); t99 * t224 + qJDD(1) * mrSges(4,1) - t130 * mrSges(4,3) + t63 * t17 + t65 * t18 + t274 * t51 + t276 * t50 + t267 * t121 + t131 * t120 + m(6) * (t120 * t133 - t121 * t34) + m(5) * (t120 * t39 + t121 * t38) + m(4) * t97 + (t245 * t120 - t121 * t268) * t243 + (t129 * g(1) + t126 * g(2)) * t184 + (t1 * t65 + t10 * t276 + t120 * t148 + t2 * t63 + t274 * t9) * m(7); m(5) * qJDD(4) + t190 * g(3) + ((-t124 * t51 + t127 * t50 + t95) * qJD(5) + m(7) * (t10 * t197 - t199 * t9 - t8) + m(6) * (t12 + t208) - t226) * t128 + (t58 + (-t124 * t50 - t127 * t51) * qJD(6) + t225 * qJD(5) + m(7) * (qJD(5) * t21 + t270) + m(6) * (-qJD(5) * t32 + t11) + t254) * t125; (t101 * t153 + t156 * t81 + t159 * t82) * qJD(6) / 0.2e1 - (t101 * t134 + t81 * t135 + t136 * t82) * qJD(1) / 0.2e1 + (t200 * t26 + t6) * t124 / 0.2e1 + t8 * t161 + t270 * mrSges(7,3) + t127 * t5 / 0.2e1 + Ifges(6,6) * t90 + Ifges(6,5) * t91 + (-t95 + t181) * t32 - t16 * t50 - t15 * t51 - t11 * mrSges(6,2) + t12 * mrSges(6,1) - pkin(5) * t13 + (-t10 * (-mrSges(7,2) * t125 - mrSges(7,3) * t205) - t9 * (mrSges(7,1) * t125 - mrSges(7,3) * t203) - t144) * qJD(1) + (-pkin(5) * t8 - t10 * t16 - t15 * t9) * m(7) + (-t163 - t132) * g(3) + t27 * t193 / 0.2e1 - t154 * t188 / 0.2e1 + (t125 * t141 - t128 * t183 + t162) * (g(1) * t80 - g(2) * t79) - (-Ifges(6,2) * t201 + t107 + t258) * t200 / 0.2e1 + (m(7) * ((-t10 * t124 - t127 * t9) * qJD(6) + t164) - t51 * t193 - t50 * t195 + t254) * pkin(8) + t250 * t201 + (t182 - t248) * t33 + Ifges(6,3) * qJDD(5) + t101 * t21 * t160 + t152 * t237 + t155 * t238 + t158 * t239 + t195 * t240 - t130 * t143 / 0.2e1; -t21 * (mrSges(7,1) * t82 + mrSges(7,2) * t81) + (Ifges(7,1) * t81 - t231) * t278 + t26 * t235 + (Ifges(7,5) * t81 - Ifges(7,6) * t82) * t277 - t9 * t50 + t10 * t51 - g(1) * (mrSges(7,1) * t30 - mrSges(7,2) * t31) - g(2) * ((-t127 * t80 + t205 * t79) * mrSges(7,1) + (t124 * t80 + t203 * t79) * mrSges(7,2)) + g(3) * t72 + (t10 * t82 + t81 * t9) * mrSges(7,3) + t186 + (-Ifges(7,2) * t82 + t27 + t75) * t279 - t246;];
tau  = t23;
