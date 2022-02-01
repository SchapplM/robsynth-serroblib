% Calculate vector of inverse dynamics joint torques for
% S5RRPRR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2022-01-20 10:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:48:03
% EndTime: 2022-01-20 10:48:12
% DurationCPUTime: 2.69s
% Computational Cost: add. (3753->363), mult. (5994->484), div. (0->0), fcn. (3531->16), ass. (0->184)
t182 = sin(qJ(4));
t234 = Ifges(5,4) * t182;
t260 = t234 / 0.2e1;
t186 = cos(qJ(4));
t233 = Ifges(5,4) * t186;
t259 = t186 * Ifges(5,2);
t177 = qJ(4) + qJ(5);
t167 = sin(t177);
t235 = mrSges(6,2) * t167;
t236 = mrSges(5,2) * t182;
t258 = t235 + t236;
t178 = qJ(1) + qJ(2);
t165 = pkin(9) + t178;
t151 = sin(t165);
t152 = cos(t165);
t253 = g(1) * t152 + g(2) * t151;
t257 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t169 = cos(t177);
t149 = t169 * mrSges(6,1);
t229 = t186 * mrSges(5,1);
t256 = -mrSges(4,1) - t149 - t229;
t179 = sin(pkin(9));
t183 = sin(qJ(2));
t232 = pkin(1) * qJD(1);
t213 = t183 * t232;
t137 = t179 * t213;
t180 = cos(pkin(9));
t187 = cos(qJ(2));
t212 = t187 * t232;
t105 = t180 * t212 - t137;
t246 = pkin(2) * t179;
t153 = pkin(7) + t246;
t237 = -pkin(8) - t153;
t205 = qJD(4) * t237;
t109 = t182 * t205;
t110 = t186 * t205;
t181 = sin(qJ(5));
t185 = cos(qJ(5));
t125 = t181 * t186 + t182 * t185;
t119 = t237 * t182;
t171 = t186 * pkin(8);
t120 = t153 * t186 + t171;
t72 = t119 * t181 + t120 * t185;
t255 = -qJD(5) * t72 + t125 * t105 - t109 * t181 + t110 * t185;
t124 = -t181 * t182 + t185 * t186;
t71 = t119 * t185 - t120 * t181;
t254 = qJD(5) * t71 - t124 * t105 + t109 * t185 + t110 * t181;
t252 = m(3) * pkin(1);
t176 = qJD(1) + qJD(2);
t99 = t125 * t176;
t250 = t99 / 0.2e1;
t98 = t124 * t176;
t249 = mrSges(6,3) * t98;
t248 = Ifges(6,4) * t99;
t247 = pkin(1) * t187;
t170 = cos(t178);
t156 = pkin(2) * t170;
t245 = pkin(2) * t180;
t244 = pkin(4) * t186;
t166 = t186 * qJD(3);
t131 = pkin(2) * t176 + t212;
t87 = t179 * t131 + t180 * t213;
t77 = pkin(7) * t176 + t87;
t67 = -t182 * t77 + t166;
t241 = t67 * mrSges(5,3);
t217 = qJD(3) * t182;
t68 = t186 * t77 + t217;
t240 = t68 * mrSges(5,3);
t239 = t99 * mrSges(6,3);
t158 = pkin(2) + t247;
t219 = t180 * t183;
t108 = pkin(1) * t219 + t179 * t158;
t102 = pkin(7) + t108;
t238 = -pkin(8) - t102;
t208 = pkin(8) * t176 + t77;
t58 = t186 * t208 + t217;
t231 = t181 * t58;
t230 = t185 * t58;
t216 = qJD(4) * t182;
t174 = qJDD(1) + qJDD(2);
t121 = -qJD(2) * t213 + qJDD(1) * t247;
t100 = pkin(2) * t174 + t121;
t218 = qJD(2) * t187;
t122 = (qJD(1) * t218 + qJDD(1) * t183) * pkin(1);
t66 = t179 * t100 + t180 * t122;
t57 = pkin(7) * t174 + t66;
t20 = qJD(4) * t166 + t182 * qJDD(3) + t186 * t57 - t216 * t77;
t228 = t186 * t20;
t227 = Ifges(5,5) * qJD(4);
t226 = Ifges(5,6) * qJD(4);
t225 = qJD(4) * t68;
t220 = t179 * t183;
t106 = (t180 * t187 - t220) * qJD(2) * pkin(1);
t224 = t106 * t182;
t223 = t106 * t186;
t222 = t176 * t182;
t221 = t176 * t186;
t215 = qJD(4) * t186;
t214 = m(4) + m(5) + m(6);
t211 = pkin(4) * t216;
t210 = t152 * pkin(3) + t151 * pkin(7) + t156;
t157 = pkin(3) + t244;
t209 = t176 * t216;
t206 = qJD(4) * t238;
t204 = t149 - t235;
t65 = t100 * t180 - t179 * t122;
t86 = t131 * t180 - t137;
t107 = -pkin(1) * t220 + t158 * t180;
t101 = -pkin(3) - t107;
t189 = -pkin(8) - pkin(7);
t203 = -t151 * t189 + t152 * t157 + t156;
t135 = -t229 + t236;
t202 = mrSges(5,1) * t182 + mrSges(5,2) * t186;
t201 = mrSges(6,1) * t167 + mrSges(6,2) * t169;
t200 = t234 + t259;
t56 = -t182 * t208 + t166;
t50 = qJD(4) * pkin(4) + t56;
t15 = t185 * t50 - t231;
t16 = t181 * t50 + t230;
t82 = t238 * t182;
t83 = t102 * t186 + t171;
t43 = -t181 * t83 + t185 * t82;
t44 = t181 * t82 + t185 * t83;
t199 = -t182 * t67 + t186 * t68;
t55 = -pkin(3) * t174 - t65;
t198 = pkin(1) * (t179 * t187 + t219);
t197 = t125 * qJD(5);
t196 = t124 * qJD(5);
t104 = qJD(2) * t198;
t21 = t186 * qJDD(3) - t182 * t57 - t225;
t168 = sin(t178);
t194 = -t170 * mrSges(3,1) + mrSges(3,2) * t168 + t257 * t151 + t256 * t152;
t173 = qJDD(4) + qJDD(5);
t175 = qJD(4) + qJD(5);
t115 = t174 * t182 + t176 * t215;
t12 = qJDD(4) * pkin(4) - pkin(8) * t115 + t21;
t114 = t174 * t186 - t209;
t13 = pkin(8) * t114 + t20;
t3 = qJD(5) * t15 + t12 * t181 + t13 * t185;
t4 = -qJD(5) * t16 + t12 * t185 - t13 * t181;
t41 = t114 * t181 + t115 * t185 + t176 * t196;
t42 = t114 * t185 - t115 * t181 - t176 * t197;
t51 = Ifges(6,2) * t98 + Ifges(6,6) * t175 + t248;
t93 = Ifges(6,4) * t98;
t52 = Ifges(6,1) * t99 + Ifges(6,5) * t175 + t93;
t69 = -t157 * t176 - t86;
t193 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t15 * t249 + t51 * t250 - t69 * (mrSges(6,1) * t99 + mrSges(6,2) * t98) + Ifges(6,3) * t173 - t99 * (Ifges(6,1) * t98 - t248) / 0.2e1 + Ifges(6,6) * t42 + Ifges(6,5) * t41 - t175 * (Ifges(6,5) * t98 - Ifges(6,6) * t99) / 0.2e1 - (-Ifges(6,2) * t99 + t52 + t93) * t98 / 0.2e1;
t192 = m(5) * (-t182 * t21 + t228 + (-t182 * t68 - t186 * t67) * qJD(4));
t191 = mrSges(3,2) * t170 + (pkin(2) * t214 + mrSges(3,1)) * t168 + (m(5) * pkin(3) + m(6) * t157 - t256 - t258) * t151 + (-m(5) * pkin(7) + m(6) * t189 + t257) * t152;
t38 = -pkin(4) * t114 + t55;
t74 = qJD(4) * t124 + t196;
t75 = -qJD(4) * t125 - t197;
t76 = -pkin(3) * t176 - t86;
t96 = t176 * t200 + t226;
t138 = Ifges(5,4) * t221;
t97 = Ifges(5,1) * t222 + t138 + t227;
t190 = t76 * t202 * qJD(4) + t175 * (Ifges(6,5) * t74 + Ifges(6,6) * t75) / 0.2e1 + t55 * t135 + t121 * mrSges(3,1) - t122 * mrSges(3,2) + t98 * (Ifges(6,4) * t74 + Ifges(6,2) * t75) / 0.2e1 + t74 * t52 / 0.2e1 + t75 * t51 / 0.2e1 + t69 * (-mrSges(6,1) * t75 + mrSges(6,2) * t74) + t65 * mrSges(4,1) - t96 * t216 / 0.2e1 + mrSges(5,3) * t228 + (Ifges(5,1) * t186 - t234) * t209 / 0.2e1 + (Ifges(6,1) * t74 + Ifges(6,4) * t75) * t250 + qJD(4) ^ 2 * (Ifges(5,5) * t186 - Ifges(5,6) * t182) / 0.2e1 + (t259 / 0.2e1 + t260 + t200 / 0.2e1) * t114 + t115 * (t182 * Ifges(5,1) + t233) + (t176 * (-Ifges(5,2) * t182 + t233) + t97) * t215 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t174 + (-t15 * t74 + t16 * t75) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t182 + Ifges(5,6) * t186) + (t38 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t41 + Ifges(6,4) * t42 + Ifges(6,5) * t173) * t125 + (-t38 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t41 + Ifges(6,2) * t42 + Ifges(6,6) * t173) * t124;
t188 = cos(qJ(1));
t184 = sin(qJ(1));
t172 = t188 * pkin(1);
t154 = -pkin(3) - t245;
t134 = -t157 - t245;
t130 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t221;
t129 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t222;
t113 = t135 * t176;
t103 = qJD(1) * t198;
t92 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t115;
t91 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t114;
t90 = t101 - t244;
t88 = t104 + t211;
t81 = mrSges(6,1) * t175 - t239;
t80 = -mrSges(6,2) * t175 + t249;
t70 = -mrSges(5,1) * t114 + mrSges(5,2) * t115;
t62 = -mrSges(6,1) * t98 + mrSges(6,2) * t99;
t47 = t186 * t206 - t224;
t46 = t182 * t206 + t223;
t35 = -mrSges(6,2) * t173 + mrSges(6,3) * t42;
t34 = mrSges(6,1) * t173 - mrSges(6,3) * t41;
t19 = t185 * t56 - t231;
t18 = -t181 * t56 - t230;
t10 = -mrSges(6,1) * t42 + mrSges(6,2) * t41;
t6 = -qJD(5) * t44 - t181 * t46 + t185 * t47;
t5 = qJD(5) * t43 + t181 * t47 + t185 * t46;
t1 = [m(4) * (-t104 * t86 + t106 * t87 + t107 * t65 + t108 * t66) + m(6) * (t15 * t6 + t16 * t5 + t3 * t44 + t38 * t90 + t4 * t43 + t69 * t88) + m(5) * (t101 * t55 + t104 * t76 + t223 * t68 - t224 * t67) + (mrSges(2,2) * t188 + (mrSges(2,1) + (m(3) + t214) * pkin(1)) * t184 + t191) * g(1) + t190 + (-t104 * t176 + t107 * t174) * mrSges(4,1) + t104 * t113 + (mrSges(2,2) * t184 - m(6) * (t172 + t203) - m(5) * (t172 + t210) - m(4) * (t156 + t172) + (-mrSges(2,1) - t252) * t188 + t258 * t152 + t194) * g(2) + t101 * t70 + t88 * t62 + t90 * t10 + t5 * t80 + t6 * t81 + (-t106 * t129 + (-qJD(4) * t130 - t92) * t102 + (-t21 - t225) * mrSges(5,3)) * t182 + t43 * t34 + t44 * t35 + (-t106 * t176 - t108 * t174 - t66) * mrSges(4,2) + ((-t174 * t183 - t176 * t218) * mrSges(3,2) + (-qJD(2) * t176 * t183 + t174 * t187) * mrSges(3,1)) * pkin(1) + (t102 * t91 + t106 * t130 + (-t102 * t129 - t241) * qJD(4)) * t186 + (t121 * t187 + t122 * t183) * t252 + t102 * t192 + Ifges(2,3) * qJDD(1); t190 + t154 * t70 + t174 * mrSges(4,1) * t245 + (t152 * t235 + t194) * g(2) + t134 * t10 + t71 * t34 + t72 * t35 + t255 * t81 + t254 * t80 + (-t174 * t246 - t66) * mrSges(4,2) + (-t113 - t62) * t103 + (-t105 * t130 + t153 * t91 + (-t129 * t153 - t241) * qJD(4)) * t186 + (mrSges(4,1) * t103 + mrSges(4,2) * t105 + (mrSges(3,1) * t183 + mrSges(3,2) * t187) * t232) * t176 + t153 * t192 + t191 * g(1) + (g(2) * mrSges(5,2) * t152 - t21 * mrSges(5,3) + t105 * t129 - t153 * t92 + (pkin(4) * t62 - t130 * t153 - t240) * qJD(4)) * t182 + (-t203 * g(2) + t134 * t38 + t3 * t72 + t4 * t71 + (-t103 + t211) * t69 + t254 * t16 + t255 * t15) * m(6) + (-g(2) * t210 - t103 * t76 - t105 * t199 + t154 * t55) * m(5) + (-t156 * g(2) + t103 * t86 - t105 * t87 + (t179 * t66 + t180 * t65) * pkin(2)) * m(4); m(4) * qJDD(3) + t124 * t34 + t125 * t35 + t182 * t91 + t186 * t92 + t74 * t80 + t75 * t81 + (-t182 * t129 + t186 * t130) * qJD(4) + m(5) * (qJD(4) * t199 + t182 * t20 + t186 * t21) + m(6) * (t124 * t4 + t125 * t3 + t15 * t75 + t16 * t74) - t214 * g(3); t193 + (t135 - t204) * g(3) + t16 * t239 + (t181 * t35 + t185 * t34 + (-g(3) * t186 + t181 * t3 + t253 * t182 + t185 * t4) * m(6) + (-t181 * t81 + t185 * t80 + (-t15 * t181 + t16 * t185) * m(6)) * qJD(5)) * pkin(4) - m(6) * (t15 * t18 + t16 * t19) + ((-t76 * mrSges(5,2) - t227 / 0.2e1 + t241 - t138 / 0.2e1 - t97 / 0.2e1) * t186 + (-t76 * mrSges(5,1) + t226 / 0.2e1 + t240 + t96 / 0.2e1 + (t260 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t186) * t176 + (-m(6) * t69 - t62) * pkin(4)) * t182) * t176 + t68 * t129 - t67 * t130 + Ifges(5,6) * t114 + Ifges(5,5) * t115 - t19 * t80 - t18 * t81 - t20 * mrSges(5,2) + t21 * mrSges(5,1) + Ifges(5,3) * qJDD(4) + t253 * (t201 + t202); t193 + (t81 + t239) * t16 - g(3) * t204 - t15 * t80 + t253 * t201;];
tau = t1;
