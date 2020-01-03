% Calculate vector of inverse dynamics joint torques for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRRR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:53:47
% EndTime: 2020-01-03 11:53:53
% DurationCPUTime: 2.40s
% Computational Cost: add. (3525->352), mult. (6055->465), div. (0->0), fcn. (3569->16), ass. (0->176)
t179 = sin(qJ(4));
t238 = Ifges(5,4) * t179;
t261 = t238 / 0.2e1;
t183 = cos(qJ(4));
t237 = Ifges(5,4) * t183;
t260 = t183 * Ifges(5,2);
t241 = mrSges(5,1) * t183;
t259 = -mrSges(4,1) - t241;
t177 = cos(pkin(9));
t154 = pkin(1) * t177 + pkin(2);
t176 = sin(pkin(9));
t247 = pkin(1) * t176;
t212 = qJD(1) * t247;
t258 = -qJD(3) * t212 + t154 * qJDD(1);
t129 = t154 * qJD(1);
t257 = qJD(3) * t129 + qJDD(1) * t247;
t256 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t255 = m(3) * pkin(1);
t178 = sin(qJ(5));
t182 = cos(qJ(5));
t114 = t178 * t183 + t179 * t182;
t186 = -pkin(8) - pkin(7);
t208 = qJD(4) * t186;
t121 = t179 * t208;
t122 = t183 * t208;
t133 = t186 * t179;
t168 = t183 * pkin(8);
t134 = pkin(7) * t183 + t168;
t86 = t133 * t178 + t134 * t182;
t180 = sin(qJ(3));
t184 = cos(qJ(3));
t88 = t129 * t184 - t180 * t212;
t254 = -qJD(5) * t86 + t114 * t88 - t121 * t178 + t122 * t182;
t113 = -t178 * t179 + t182 * t183;
t85 = t133 * t182 - t134 * t178;
t253 = qJD(5) * t85 - t113 * t88 + t121 * t182 + t122 * t178;
t103 = t154 * t184 - t180 * t247;
t104 = t180 * t154 + t184 * t247;
t252 = m(5) * pkin(7);
t250 = m(3) + m(4);
t173 = qJD(1) + qJD(3);
t100 = t114 * t173;
t249 = t100 / 0.2e1;
t99 = t113 * t173;
t248 = mrSges(6,3) * t99;
t246 = pkin(4) * t183;
t174 = qJ(1) + pkin(9);
t163 = qJ(3) + t174;
t152 = sin(t163);
t245 = g(2) * t152;
t164 = t183 * qJD(2);
t89 = t129 * t180 + t184 * t212;
t77 = pkin(7) * t173 + t89;
t68 = -t179 * t77 + t164;
t244 = t68 * mrSges(5,3);
t216 = qJD(2) * t179;
t69 = t183 * t77 + t216;
t243 = t69 * mrSges(5,3);
t102 = pkin(7) + t104;
t242 = -pkin(8) - t102;
t240 = mrSges(5,2) * t179;
t175 = qJ(4) + qJ(5);
t165 = sin(t175);
t239 = mrSges(6,2) * t165;
t236 = Ifges(6,4) * t100;
t235 = t100 * mrSges(6,3);
t166 = cos(t175);
t151 = t166 * mrSges(6,1);
t214 = qJD(4) * t179;
t171 = qJDD(1) + qJDD(3);
t57 = t180 * t258 + t184 * t257;
t51 = pkin(7) * t171 + t57;
t17 = qJD(4) * t164 + qJDD(2) * t179 + t183 * t51 - t214 * t77;
t234 = t17 * t183;
t206 = pkin(8) * t173 + t77;
t63 = t183 * t206 + t216;
t233 = t178 * t63;
t95 = t103 * qJD(3);
t232 = t179 * t95;
t18 = -qJD(4) * t69 + qJDD(2) * t183 - t179 * t51;
t231 = t18 * t179;
t230 = t182 * t63;
t229 = t183 * t95;
t228 = Ifges(5,5) * qJD(4);
t227 = Ifges(5,6) * qJD(4);
t153 = cos(t163);
t226 = t153 * t165;
t225 = t153 * t166;
t223 = t173 * t179;
t222 = t173 * t183;
t156 = pkin(3) + t246;
t221 = t152 * t156 + t153 * t186;
t220 = mrSges(6,1) * t226 + mrSges(6,2) * t225;
t219 = pkin(3) * t153 + pkin(7) * t152;
t159 = sin(t174);
t181 = sin(qJ(1));
t218 = pkin(1) * t181 + pkin(2) * t159;
t160 = cos(t174);
t185 = cos(qJ(1));
t217 = t185 * pkin(1) + pkin(2) * t160;
t213 = qJD(4) * t183;
t211 = pkin(4) * t214;
t210 = -mrSges(2,1) - t255;
t207 = t173 * t214;
t204 = qJD(4) * t242;
t203 = t151 - t239;
t202 = -t152 * t186 + t153 * t156;
t101 = -pkin(3) - t103;
t130 = t240 - t241;
t200 = mrSges(5,1) * t179 + mrSges(5,2) * t183;
t199 = -mrSges(6,1) * t165 - mrSges(6,2) * t166;
t198 = t239 + t240;
t197 = t238 + t260;
t62 = -t179 * t206 + t164;
t59 = qJD(4) * pkin(4) + t62;
t19 = t182 * t59 - t233;
t20 = t178 * t59 + t230;
t80 = t242 * t179;
t81 = t102 * t183 + t168;
t42 = -t178 * t81 + t182 * t80;
t43 = t178 * t80 + t182 * t81;
t196 = -t179 * t68 + t183 * t69;
t123 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t223;
t124 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t222;
t195 = -t179 * t123 + t183 * t124;
t58 = -t180 * t257 + t184 * t258;
t194 = t114 * qJD(5);
t193 = t113 * qJD(5);
t192 = -mrSges(6,1) * t225 + t152 * t256 + t153 * t259;
t96 = t104 * qJD(3);
t52 = -pkin(3) * t171 - t58;
t191 = -t231 + (-t179 * t69 - t183 * t68) * qJD(4);
t190 = (t252 - t256) * t153 + (-t151 + t259) * t152;
t170 = qJDD(4) + qJDD(5);
t172 = qJD(4) + qJD(5);
t107 = t171 * t179 + t173 * t213;
t11 = qJDD(4) * pkin(4) - pkin(8) * t107 + t18;
t106 = t171 * t183 - t207;
t13 = pkin(8) * t106 + t17;
t3 = qJD(5) * t19 + t11 * t178 + t13 * t182;
t4 = -qJD(5) * t20 + t11 * t182 - t13 * t178;
t40 = t106 * t178 + t107 * t182 + t173 * t193;
t41 = t106 * t182 - t107 * t178 - t173 * t194;
t55 = Ifges(6,2) * t99 + Ifges(6,6) * t172 + t236;
t92 = Ifges(6,4) * t99;
t56 = Ifges(6,1) * t100 + Ifges(6,5) * t172 + t92;
t71 = -t156 * t173 - t88;
t189 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t19 * t248 + t55 * t249 - t71 * (mrSges(6,1) * t100 + mrSges(6,2) * t99) + Ifges(6,3) * t170 - t100 * (Ifges(6,1) * t99 - t236) / 0.2e1 + Ifges(6,6) * t41 + Ifges(6,5) * t40 - t172 * (Ifges(6,5) * t99 - Ifges(6,6) * t100) / 0.2e1 - (-Ifges(6,2) * t100 + t56 + t92) * t99 / 0.2e1;
t31 = -pkin(4) * t106 + t52;
t72 = qJD(4) * t113 + t193;
t73 = -qJD(4) * t114 - t194;
t76 = -pkin(3) * t173 - t88;
t97 = t173 * t197 + t227;
t138 = Ifges(5,4) * t222;
t98 = Ifges(5,1) * t223 + t138 + t228;
t188 = t76 * t200 * qJD(4) + (Ifges(6,1) * t72 + Ifges(6,4) * t73) * t249 + mrSges(5,3) * t234 + t172 * (Ifges(6,5) * t72 + Ifges(6,6) * t73) / 0.2e1 + Ifges(4,3) * t171 + t52 * t130 + t99 * (Ifges(6,4) * t72 + Ifges(6,2) * t73) / 0.2e1 + t72 * t56 / 0.2e1 + t73 * t55 / 0.2e1 + t71 * (-mrSges(6,1) * t73 + mrSges(6,2) * t72) + t58 * mrSges(4,1) + (Ifges(5,1) * t183 - t238) * t207 / 0.2e1 - t97 * t214 / 0.2e1 + qJD(4) ^ 2 * (Ifges(5,5) * t183 - Ifges(5,6) * t179) / 0.2e1 + (t197 / 0.2e1 + t260 / 0.2e1 + t261) * t106 + t107 * (t179 * Ifges(5,1) + t237) + (t173 * (-Ifges(5,2) * t179 + t237) + t98) * t213 / 0.2e1 + (-t19 * t72 + t20 * t73) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t179 + Ifges(5,6) * t183) + (mrSges(6,2) * t31 - mrSges(6,3) * t4 + Ifges(6,1) * t40 + Ifges(6,4) * t41 + Ifges(6,5) * t170) * t114 + (-mrSges(6,1) * t31 + mrSges(6,3) * t3 + Ifges(6,4) * t40 + Ifges(6,2) * t41 + Ifges(6,6) * t170) * t113;
t145 = t152 * pkin(3);
t105 = t130 * t173;
t91 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t107;
t90 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t106;
t87 = t101 - t246;
t82 = t96 + t211;
t79 = mrSges(6,1) * t172 - t235;
t78 = -mrSges(6,2) * t172 + t248;
t67 = -mrSges(5,1) * t106 + mrSges(5,2) * t107;
t65 = -mrSges(6,1) * t99 + mrSges(6,2) * t100;
t50 = t183 * t204 - t232;
t49 = t179 * t204 + t229;
t34 = -mrSges(6,2) * t170 + mrSges(6,3) * t41;
t33 = mrSges(6,1) * t170 - mrSges(6,3) * t40;
t22 = t182 * t62 - t233;
t21 = -t178 * t62 - t230;
t10 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t6 = -qJD(5) * t43 - t178 * t49 + t182 * t50;
t5 = qJD(5) * t42 + t178 * t50 + t182 * t49;
t1 = [t195 * t95 + t188 + t96 * t105 + t101 * t67 + m(5) * (t101 * t52 + t69 * t229 - t68 * t232 + t76 * t96) + (mrSges(2,2) * t181 - mrSges(3,1) * t160 + mrSges(3,2) * t159 - m(6) * (t202 + t217) - m(5) * (t217 + t219) - m(4) * t217 + t210 * t185 + t198 * t153 + t192) * g(2) + t82 * t65 + t87 * t10 + t5 * t78 + t6 * t79 + t42 * t33 + t43 * t34 + (-mrSges(2,2) * t185 - m(5) * (t145 + t218) - mrSges(3,1) * t159 - mrSges(3,2) * t160 - m(6) * (t218 + t221) - m(4) * t218 + t210 * t181 + t198 * t152 + t190) * g(3) + t191 * mrSges(5,3) + (t103 * t171 - t173 * t96) * mrSges(4,1) + (-t104 * t171 - t173 * t95 - t57) * mrSges(4,2) + m(4) * (t103 * t58 + t104 * t57 - t88 * t96 + t89 * t95) + m(6) * (t19 * t6 + t20 * t5 + t3 * t43 + t31 * t87 + t4 * t42 + t71 * t82) + (m(5) * (-t213 * t68 - t214 * t69 - t231 + t234) + t183 * t90 - t179 * t91 - t124 * t214 - t123 * t213) * t102 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * mrSges(3,1) * t177 - 0.2e1 * mrSges(3,2) * t176 + (t176 ^ 2 + t177 ^ 2) * t255) * pkin(1)) * qJDD(1); t113 * t33 + t114 * t34 + t179 * t90 + t183 * t91 + t72 * t78 + t73 * t79 + t250 * qJDD(2) + t195 * qJD(4) + m(5) * (qJD(4) * t196 + t17 * t179 + t18 * t183) + m(6) * (t113 * t4 + t114 * t3 + t19 * t73 + t20 * t72) + (-m(5) - m(6) - t250) * g(1); (t191 + t234) * t252 + (t152 * t239 + t190) * g(3) + t188 + (t173 * t88 - t57) * mrSges(4,2) + (mrSges(6,2) * t226 + t192) * g(2) - t156 * t10 + t85 * t33 + t86 * t34 - pkin(3) * t67 + (t173 * mrSges(4,1) - t105 - t65) * t89 + t254 * t79 + t253 * t78 + (-qJD(4) * t244 - t88 * t124 + (-qJD(4) * t123 + t90) * pkin(7)) * t183 + (-t18 * mrSges(5,3) - pkin(7) * t91 + t88 * t123 + (g(2) * t153 + g(3) * t152) * mrSges(5,2) + (pkin(4) * t65 - pkin(7) * t124 - t243) * qJD(4)) * t179 + (-t202 * g(2) - t221 * g(3) - t156 * t31 + t3 * t86 + t4 * t85 + (t211 - t89) * t71 + t253 * t20 + t254 * t19) * m(6) + (-pkin(3) * t52 - g(2) * t219 - g(3) * t145 - t196 * t88 - t76 * t89) * m(5); t189 + t20 * t235 + t69 * t123 - t68 * t124 + Ifges(5,6) * t106 + Ifges(5,5) * t107 - t22 * t78 - t21 * t79 - t17 * mrSges(5,2) + t18 * mrSges(5,1) - m(6) * (t19 * t21 + t20 * t22) + (-t199 + t200) * t245 + (t130 - t203) * g(1) + Ifges(5,3) * qJDD(4) + (t178 * t34 + t182 * t33 + (t178 * t3 + t182 * t4 - g(1) * t183 + (-g(3) * t153 + t245) * t179) * m(6) + (-t178 * t79 + t182 * t78 + (-t178 * t19 + t182 * t20) * m(6)) * qJD(5)) * pkin(4) + ((-t138 / 0.2e1 - t98 / 0.2e1 - t76 * mrSges(5,2) - t228 / 0.2e1 + t244) * t183 + (t97 / 0.2e1 - t76 * mrSges(5,1) + t227 / 0.2e1 + t243 + (t261 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t183) * t173 + (-m(6) * t71 - t65) * pkin(4)) * t179) * t173 + (-t153 * t200 - t220) * g(3); t189 - t199 * t245 + (t79 + t235) * t20 - g(1) * t203 - g(3) * t220 - t19 * t78;];
tau = t1;
