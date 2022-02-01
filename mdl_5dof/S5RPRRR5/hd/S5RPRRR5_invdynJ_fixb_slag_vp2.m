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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:48:35
% EndTime: 2022-01-20 09:48:43
% DurationCPUTime: 2.43s
% Computational Cost: add. (3525->344), mult. (6055->453), div. (0->0), fcn. (3569->16), ass. (0->171)
t173 = sin(qJ(4));
t230 = Ifges(5,4) * t173;
t254 = t230 / 0.2e1;
t177 = cos(qJ(4));
t229 = Ifges(5,4) * t177;
t253 = t177 * Ifges(5,2);
t169 = qJ(4) + qJ(5);
t161 = cos(t169);
t146 = t161 * mrSges(6,1);
t252 = -mrSges(4,1) - t146;
t171 = cos(pkin(9));
t149 = pkin(1) * t171 + pkin(2);
t170 = sin(pkin(9));
t239 = pkin(1) * t170;
t207 = qJD(1) * t239;
t251 = -qJD(3) * t207 + t149 * qJDD(1);
t126 = t149 * qJD(1);
t250 = qJD(3) * t126 + qJDD(1) * t239;
t160 = sin(t169);
t231 = mrSges(6,2) * t160;
t232 = mrSges(5,2) * t173;
t249 = t231 + t232;
t168 = qJ(1) + pkin(9);
t158 = qJ(3) + t168;
t147 = sin(t158);
t148 = cos(t158);
t244 = g(1) * t148 + g(2) * t147;
t248 = mrSges(4,2) - mrSges(5,3) - mrSges(6,3);
t247 = m(3) * pkin(1);
t172 = sin(qJ(5));
t176 = cos(qJ(5));
t113 = -t172 * t173 + t176 * t177;
t180 = -pkin(8) - pkin(7);
t203 = qJD(4) * t180;
t118 = t173 * t203;
t119 = t177 * t203;
t129 = t180 * t173;
t162 = t177 * pkin(8);
t130 = pkin(7) * t177 + t162;
t85 = t129 * t176 - t130 * t172;
t174 = sin(qJ(3));
t178 = cos(qJ(3));
t88 = t126 * t178 - t174 * t207;
t246 = qJD(5) * t85 - t113 * t88 + t118 * t176 + t119 * t172;
t114 = t172 * t177 + t173 * t176;
t86 = t129 * t172 + t130 * t176;
t245 = -qJD(5) * t86 + t114 * t88 - t118 * t172 + t119 * t176;
t103 = t149 * t178 - t174 * t239;
t104 = t174 * t149 + t178 * t239;
t243 = m(5) * pkin(3);
t167 = qJD(1) + qJD(3);
t100 = t114 * t167;
t241 = t100 / 0.2e1;
t99 = t113 * t167;
t240 = mrSges(6,3) * t99;
t238 = pkin(4) * t177;
t159 = t177 * qJD(2);
t89 = t126 * t174 + t178 * t207;
t77 = pkin(7) * t167 + t89;
t68 = -t173 * t77 + t159;
t235 = t68 * mrSges(5,3);
t212 = qJD(2) * t173;
t69 = t177 * t77 + t212;
t234 = t69 * mrSges(5,3);
t102 = pkin(7) + t104;
t233 = -pkin(8) - t102;
t228 = Ifges(6,4) * t100;
t227 = t100 * mrSges(6,3);
t210 = qJD(4) * t173;
t165 = qJDD(1) + qJDD(3);
t57 = t251 * t174 + t250 * t178;
t51 = pkin(7) * t165 + t57;
t17 = qJD(4) * t159 + t173 * qJDD(2) + t177 * t51 - t210 * t77;
t226 = t17 * t177;
t201 = pkin(8) * t167 + t77;
t63 = t177 * t201 + t212;
t225 = t172 * t63;
t95 = t103 * qJD(3);
t224 = t173 * t95;
t223 = t176 * t63;
t222 = t177 * mrSges(5,1);
t221 = t177 * t95;
t18 = -qJD(4) * t69 + t177 * qJDD(2) - t173 * t51;
t220 = t18 * t173;
t219 = Ifges(5,5) * qJD(4);
t218 = Ifges(5,6) * qJD(4);
t216 = t167 * t173;
t215 = t167 * t177;
t214 = t148 * pkin(3) + t147 * pkin(7);
t155 = cos(t168);
t179 = cos(qJ(1));
t213 = t179 * pkin(1) + pkin(2) * t155;
t209 = qJD(4) * t177;
t208 = m(4) + m(5) + m(6);
t206 = pkin(4) * t210;
t205 = m(3) + t208;
t151 = pkin(3) + t238;
t202 = t167 * t210;
t199 = qJD(4) * t233;
t198 = t146 - t231;
t197 = -t147 * t180 + t148 * t151;
t101 = -pkin(3) - t103;
t127 = -t222 + t232;
t195 = mrSges(5,1) * t173 + mrSges(5,2) * t177;
t194 = mrSges(6,1) * t160 + mrSges(6,2) * t161;
t193 = t230 + t253;
t62 = -t173 * t201 + t159;
t59 = qJD(4) * pkin(4) + t62;
t19 = t176 * t59 - t225;
t20 = t172 * t59 + t223;
t80 = t233 * t173;
t81 = t102 * t177 + t162;
t42 = -t172 * t81 + t176 * t80;
t43 = t172 * t80 + t176 * t81;
t192 = -t173 * t68 + t177 * t69;
t120 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t216;
t121 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t215;
t191 = -t173 * t120 + t177 * t121;
t58 = -t250 * t174 + t251 * t178;
t190 = t114 * qJD(5);
t189 = t113 * qJD(5);
t188 = (-t222 + t252) * t148 + t248 * t147;
t96 = t104 * qJD(3);
t187 = m(6) * t151 + t243 - t252;
t52 = -pkin(3) * t165 - t58;
t186 = -t220 + (-t173 * t69 - t177 * t68) * qJD(4);
t184 = -t249 * t147 + (-m(5) * pkin(7) + m(6) * t180 + t248) * t148;
t164 = qJDD(4) + qJDD(5);
t166 = qJD(4) + qJD(5);
t107 = t165 * t173 + t167 * t209;
t11 = qJDD(4) * pkin(4) - pkin(8) * t107 + t18;
t106 = t165 * t177 - t202;
t13 = pkin(8) * t106 + t17;
t3 = qJD(5) * t19 + t11 * t172 + t13 * t176;
t4 = -qJD(5) * t20 + t11 * t176 - t13 * t172;
t40 = t106 * t172 + t107 * t176 + t167 * t189;
t41 = t106 * t176 - t107 * t172 - t167 * t190;
t55 = Ifges(6,2) * t99 + Ifges(6,6) * t166 + t228;
t92 = Ifges(6,4) * t99;
t56 = Ifges(6,1) * t100 + Ifges(6,5) * t166 + t92;
t71 = -t151 * t167 - t88;
t183 = t4 * mrSges(6,1) - t3 * mrSges(6,2) + t19 * t240 + t55 * t241 - t71 * (mrSges(6,1) * t100 + mrSges(6,2) * t99) + Ifges(6,3) * t164 - t100 * (Ifges(6,1) * t99 - t228) / 0.2e1 + Ifges(6,6) * t41 + Ifges(6,5) * t40 - t166 * (Ifges(6,5) * t99 - Ifges(6,6) * t100) / 0.2e1 - (-Ifges(6,2) * t100 + t56 + t92) * t99 / 0.2e1;
t31 = -pkin(4) * t106 + t52;
t72 = qJD(4) * t113 + t189;
t73 = -qJD(4) * t114 - t190;
t76 = -pkin(3) * t167 - t88;
t97 = t167 * t193 + t218;
t134 = Ifges(5,4) * t215;
t98 = Ifges(5,1) * t216 + t134 + t219;
t182 = -t97 * t210 / 0.2e1 + (Ifges(6,1) * t72 + Ifges(6,4) * t73) * t241 + mrSges(5,3) * t226 + Ifges(4,3) * t165 + t166 * (Ifges(6,5) * t72 + Ifges(6,6) * t73) / 0.2e1 + t52 * t127 + t99 * (Ifges(6,4) * t72 + Ifges(6,2) * t73) / 0.2e1 + t72 * t56 / 0.2e1 + t73 * t55 / 0.2e1 + t71 * (-mrSges(6,1) * t73 + mrSges(6,2) * t72) + t58 * mrSges(4,1) + t76 * t195 * qJD(4) + (Ifges(5,1) * t177 - t230) * t202 / 0.2e1 + qJD(4) ^ 2 * (Ifges(5,5) * t177 - Ifges(5,6) * t173) / 0.2e1 + (t193 / 0.2e1 + t253 / 0.2e1 + t254) * t106 + t107 * (t173 * Ifges(5,1) + t229) + (t98 + t167 * (-Ifges(5,2) * t173 + t229)) * t209 / 0.2e1 + (-t19 * t72 + t20 * t73) * mrSges(6,3) + qJDD(4) * (Ifges(5,5) * t173 + Ifges(5,6) * t177) + (t31 * mrSges(6,2) - t4 * mrSges(6,3) + Ifges(6,1) * t40 + Ifges(6,4) * t41 + Ifges(6,5) * t164) * t114 + (-t31 * mrSges(6,1) + t3 * mrSges(6,3) + Ifges(6,4) * t40 + Ifges(6,2) * t41 + Ifges(6,6) * t164) * t113;
t175 = sin(qJ(1));
t154 = sin(t168);
t105 = t127 * t167;
t91 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t107;
t90 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t106;
t87 = t101 - t238;
t82 = t96 + t206;
t79 = mrSges(6,1) * t166 - t227;
t78 = -mrSges(6,2) * t166 + t240;
t67 = -mrSges(5,1) * t106 + mrSges(5,2) * t107;
t65 = -mrSges(6,1) * t99 + mrSges(6,2) * t100;
t50 = t177 * t199 - t224;
t49 = t173 * t199 + t221;
t34 = -mrSges(6,2) * t164 + mrSges(6,3) * t41;
t33 = mrSges(6,1) * t164 - mrSges(6,3) * t40;
t22 = t176 * t62 - t225;
t21 = -t172 * t62 - t223;
t10 = -mrSges(6,1) * t41 + mrSges(6,2) * t40;
t6 = -qJD(5) * t43 - t172 * t49 + t176 * t50;
t5 = qJD(5) * t42 + t172 * t50 + t176 * t49;
t1 = [(m(5) * (-t209 * t68 - t210 * t69 - t220 + t226) + t177 * t90 - t173 * t91 - t120 * t209 - t121 * t210) * t102 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t171 * mrSges(3,1) - 0.2e1 * t170 * mrSges(3,2) + (t170 ^ 2 + t171 ^ 2) * t247) * pkin(1)) * qJDD(1) + t186 * mrSges(5,3) + t191 * t95 + m(4) * (t103 * t58 + t104 * t57 - t88 * t96 + t89 * t95) + m(6) * (t19 * t6 + t20 * t5 + t3 * t43 + t31 * t87 + t4 * t42 + t71 * t82) + t101 * t67 + t96 * t105 + t5 * t78 + t6 * t79 + t82 * t65 + t87 * t10 + t42 * t33 + t43 * t34 + m(5) * (t101 * t52 + t69 * t221 - t68 * t224 + t76 * t96) + (-m(6) * (t197 + t213) + mrSges(2,2) * t175 - mrSges(3,1) * t155 + mrSges(3,2) * t154 - m(5) * (t213 + t214) - m(4) * t213 + (-mrSges(2,1) - t247) * t179 + t249 * t148 + t188) * g(2) + (t103 * t165 - t167 * t96) * mrSges(4,1) + t182 + (-t104 * t165 - t167 * t95 - t57) * mrSges(4,2) + (mrSges(2,2) * t179 + mrSges(3,2) * t155 + (pkin(2) * t208 + mrSges(3,1)) * t154 + (pkin(1) * t205 + mrSges(2,1)) * t175 + (t187 + t222) * t147 + t184) * g(1); t113 * t33 + t114 * t34 + t173 * t90 + t177 * t91 + t72 * t78 + t73 * t79 + (m(3) + m(4)) * qJDD(2) + t191 * qJD(4) + m(5) * (qJD(4) * t192 + t17 * t173 + t177 * t18) + m(6) * (t113 * t4 + t114 * t3 + t19 * t73 + t20 * t72) - t205 * g(3); (t148 * t231 + t188) * g(2) + (t167 * mrSges(4,1) - t105 - t65) * t89 - t151 * t10 + t85 * t33 + t86 * t34 - pkin(3) * t67 + t182 - t52 * t243 + t245 * t79 + t246 * t78 + (g(2) * mrSges(5,2) * t148 - t18 * mrSges(5,3) - pkin(7) * t91 + t88 * t120 + (pkin(4) * t65 - pkin(7) * t121 - t234) * qJD(4)) * t173 + (g(1) * mrSges(5,1) * t147 - qJD(4) * t235 - t88 * t121 + (-qJD(4) * t120 + t90) * pkin(7)) * t177 + (t147 * t187 + t184) * g(1) + (t167 * t88 - t57) * mrSges(4,2) + (-t197 * g(2) - t151 * t31 + t3 * t86 + t4 * t85 + (t206 - t89) * t71 + t246 * t20 + t245 * t19) * m(6) + (-t214 * g(2) - t192 * t88 - t76 * t89 + (t186 + t226) * pkin(7)) * m(5); ((-t134 / 0.2e1 - t98 / 0.2e1 - t76 * mrSges(5,2) - t219 / 0.2e1 + t235) * t177 + (t97 / 0.2e1 - t76 * mrSges(5,1) + t218 / 0.2e1 + t234 + (t254 + (Ifges(5,2) / 0.2e1 - Ifges(5,1) / 0.2e1) * t177) * t167 + (-m(6) * t71 - t65) * pkin(4)) * t173) * t167 + (t172 * t34 + t176 * t33 + (-g(3) * t177 + t172 * t3 + t244 * t173 + t176 * t4) * m(6) + (-t172 * t79 + t176 * t78 + (-t172 * t19 + t176 * t20) * m(6)) * qJD(5)) * pkin(4) + t69 * t120 - t68 * t121 + Ifges(5,6) * t106 + Ifges(5,5) * t107 - t22 * t78 - t21 * t79 + (t127 - t198) * g(3) - t17 * mrSges(5,2) + t18 * mrSges(5,1) + t183 + t20 * t227 - m(6) * (t19 * t21 + t20 * t22) + Ifges(5,3) * qJDD(4) + t244 * (t194 + t195); (t79 + t227) * t20 - g(3) * t198 - t19 * t78 + t183 + t244 * t194;];
tau = t1;
