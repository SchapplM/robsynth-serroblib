% Calculate vector of inverse dynamics joint torques for
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
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
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S6RPPPRR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR1_invdynJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:29:51
% EndTime: 2019-03-09 01:29:59
% DurationCPUTime: 5.02s
% Computational Cost: add. (2346->380), mult. (4076->517), div. (0->0), fcn. (2047->10), ass. (0->181)
t92 = sin(qJ(5));
t237 = -t92 / 0.2e1;
t95 = cos(qJ(5));
t174 = qJD(1) * t95;
t151 = mrSges(6,3) * t174;
t94 = cos(qJ(6));
t169 = qJD(5) * t94;
t91 = sin(qJ(6));
t58 = -t174 * t91 + t169;
t171 = qJD(5) * t91;
t59 = t174 * t94 + t171;
t177 = -qJD(5) * mrSges(6,1) - mrSges(7,1) * t58 + mrSges(7,2) * t59 + t151;
t89 = sin(pkin(9));
t76 = pkin(1) * t89 + qJ(3);
t66 = t76 * qJD(1);
t65 = qJD(4) + t66;
t53 = -qJD(1) * pkin(7) + t65;
t34 = -qJD(2) * t92 + t53 * t95;
t29 = -qJD(5) * pkin(5) - t34;
t236 = m(7) * t29 + t177;
t235 = -t58 / 0.2e1;
t234 = -t59 / 0.2e1;
t175 = qJD(1) * t92;
t72 = qJD(6) + t175;
t233 = -t72 / 0.2e1;
t195 = Ifges(6,4) * t95;
t231 = Ifges(6,2) * t237 + t195 / 0.2e1;
t205 = t95 / 0.2e1;
t230 = -m(5) - m(4);
t229 = m(6) + m(5);
t228 = -m(7) - m(6);
t160 = qJD(1) * qJD(5);
t62 = qJDD(1) * t95 - t160 * t92;
t24 = qJD(6) * t58 + qJDD(5) * t91 + t62 * t94;
t25 = -qJD(6) * t59 + qJDD(5) * t94 - t62 * t91;
t5 = -mrSges(7,1) * t25 + mrSges(7,2) * t24;
t197 = -qJDD(5) * mrSges(6,1) + mrSges(6,3) * t62 + t5;
t227 = mrSges(4,2) - mrSges(5,3);
t226 = -mrSges(4,3) - mrSges(5,2);
t168 = qJD(5) * t95;
t170 = qJD(5) * t92;
t176 = pkin(1) * qJDD(1);
t79 = t89 * t176;
t87 = qJD(1) * qJD(3);
t55 = -qJDD(1) * qJ(3) - t79 - t87;
t47 = qJDD(4) - t55;
t42 = -qJDD(1) * pkin(7) + t47;
t14 = -qJD(2) * t170 + qJDD(2) * t95 + t168 * t53 + t42 * t92;
t225 = -qJD(5) * t34 + t14;
t224 = qJDD(1) * t76 + t87;
t35 = qJD(2) * t95 + t53 * t92;
t163 = t35 * qJD(5);
t15 = -qJDD(2) * t92 + t42 * t95 - t163;
t223 = -t163 - t15;
t196 = Ifges(6,4) * t92;
t129 = Ifges(6,1) * t95 - t196;
t52 = Ifges(7,4) * t58;
t21 = Ifges(7,1) * t59 + Ifges(7,5) * t72 + t52;
t222 = Ifges(6,5) * qJD(5) + qJD(1) * t129 + t94 * t21;
t70 = -pkin(7) + t76;
t221 = qJD(3) * t92 + t168 * t70;
t166 = qJD(6) * t94;
t167 = qJD(6) * t91;
t30 = qJD(5) * pkin(8) + t35;
t139 = pkin(5) * t92 - pkin(8) * t95;
t178 = pkin(2) + qJ(4);
t109 = -t139 - t178;
t90 = cos(pkin(9));
t202 = pkin(1) * t90;
t46 = -t109 + t202;
t31 = qJD(1) * t46 - qJD(3);
t8 = -t30 * t91 + t31 * t94;
t9 = t30 * t94 + t31 * t91;
t219 = -t166 * t8 - t167 * t9;
t63 = -qJDD(1) * t92 - t160 * t95;
t56 = qJDD(6) - t63;
t12 = mrSges(7,1) * t56 - mrSges(7,3) * t24;
t13 = -mrSges(7,2) * t56 + mrSges(7,3) * t25;
t218 = -t91 * t12 + t94 * t13;
t85 = qJ(1) + pkin(9);
t81 = sin(t85);
t82 = cos(t85);
t217 = -g(1) * t82 - g(2) * t81;
t216 = Ifges(6,6) * qJD(5) / 0.2e1 + qJD(1) * t231 + Ifges(7,5) * t234 + Ifges(7,6) * t235 + Ifges(7,3) * t233;
t215 = mrSges(3,2) + mrSges(6,3) + t226;
t130 = mrSges(7,1) * t91 + mrSges(7,2) * t94;
t194 = Ifges(7,4) * t59;
t20 = Ifges(7,2) * t58 + Ifges(7,6) * t72 + t194;
t206 = -t91 / 0.2e1;
t214 = t130 * t29 + t20 * t206;
t132 = mrSges(6,1) * t92 + mrSges(6,2) * t95;
t181 = t95 * mrSges(7,3);
t213 = mrSges(3,1) + t132 - t181 - t227;
t11 = -qJDD(5) * pkin(5) - t15;
t32 = -mrSges(7,2) * t72 + mrSges(7,3) * t58;
t33 = mrSges(7,1) * t72 - mrSges(7,3) * t59;
t120 = -t32 * t94 + t33 * t91;
t155 = mrSges(6,3) * t175;
t67 = -qJD(5) * mrSges(6,2) - t155;
t114 = -t120 + t67;
t212 = -m(6) * t223 - m(7) * (-t169 * t9 + t171 * t8 + t11) + qJD(5) * t114 - t197;
t97 = qJD(1) ^ 2;
t211 = t24 / 0.2e1;
t210 = t25 / 0.2e1;
t209 = t56 / 0.2e1;
t207 = t59 / 0.2e1;
t93 = sin(qJ(1));
t201 = pkin(1) * t93;
t10 = qJDD(5) * pkin(8) + t14;
t161 = qJD(1) * qJD(4);
t78 = -pkin(2) - t202;
t71 = qJ(4) - t78;
t43 = qJDD(1) * t71 - qJDD(3) + t161;
t16 = -pkin(5) * t63 - pkin(8) * t62 + t43;
t2 = -qJD(6) * t9 - t10 * t91 + t16 * t94;
t198 = t2 * t91;
t96 = cos(qJ(1));
t84 = t96 * pkin(1);
t193 = Ifges(7,4) * t91;
t192 = Ifges(7,4) * t94;
t191 = t70 * t92;
t188 = t91 * t92;
t186 = t92 * t94;
t185 = t94 * mrSges(7,3);
t165 = qJD(6) * t95;
t159 = Ifges(7,5) * t24 + Ifges(7,6) * t25 + Ifges(7,3) * t56;
t158 = pkin(2) * t82 + qJ(3) * t81 + t84;
t152 = t91 * t170;
t149 = m(4) + m(7) + t229;
t148 = m(7) * pkin(8) + mrSges(7,3);
t147 = qJ(3) * t82 - t201;
t143 = -t160 / 0.2e1;
t142 = qJ(4) * t82 + t158;
t140 = pkin(5) * t95 + pkin(8) * t92;
t137 = qJD(5) * t140 - qJD(6) * t191 + qJD(4);
t1 = qJD(6) * t8 + t10 * t94 + t16 * t91;
t136 = t1 * t94 - t198;
t135 = t8 * t94 + t9 * t91;
t134 = -t63 * mrSges(6,1) + t62 * mrSges(6,2);
t133 = mrSges(6,1) * t95 - mrSges(6,2) * t92;
t131 = -mrSges(7,1) * t94 + mrSges(7,2) * t91;
t128 = Ifges(7,1) * t94 - t193;
t127 = Ifges(7,1) * t91 + t192;
t125 = -Ifges(7,2) * t91 + t192;
t124 = Ifges(7,2) * t94 + t193;
t123 = -Ifges(6,5) * t92 - Ifges(6,6) * t95;
t122 = Ifges(7,5) * t94 - Ifges(7,6) * t91;
t121 = Ifges(7,5) * t91 + Ifges(7,6) * t94;
t119 = -t91 * t32 - t94 * t33;
t118 = t34 * t95 + t35 * t92;
t54 = qJD(1) * t71 - qJD(3);
t117 = qJD(4) * t54 + t43 * t71;
t112 = t54 * t133;
t111 = t92 * (-Ifges(6,2) * t95 - t196);
t110 = t95 * (-Ifges(6,1) * t92 - t195);
t108 = m(7) * pkin(5) - t131;
t107 = t165 * t91 + t169 * t92;
t106 = -t165 * t94 + t152;
t104 = qJD(6) * t46 + t221;
t103 = Ifges(7,5) * t95 - t128 * t92;
t102 = Ifges(7,6) * t95 - t125 * t92;
t101 = Ifges(7,3) * t95 - t122 * t92;
t45 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t63;
t99 = t45 + t119 * qJD(6) + t177 * qJD(5) + m(6) * t225 + m(7) * (qJD(5) * t29 + t136 + t219) + t218;
t64 = qJDD(1) * t78 + qJDD(3);
t61 = t140 * qJD(1);
t60 = t132 * qJD(1);
t51 = t130 * t95;
t40 = t186 * t82 - t81 * t91;
t39 = -t188 * t82 - t81 * t94;
t38 = -t186 * t81 - t82 * t91;
t37 = t188 * t81 - t82 * t94;
t27 = t186 * t70 + t46 * t91;
t26 = -t188 * t70 + t46 * t94;
t18 = t34 * t94 + t61 * t91;
t17 = -t34 * t91 + t61 * t94;
t7 = -t104 * t91 + t137 * t94;
t6 = t104 * t94 + t137 * t91;
t4 = Ifges(7,1) * t24 + Ifges(7,4) * t25 + Ifges(7,5) * t56;
t3 = Ifges(7,4) * t24 + Ifges(7,2) * t25 + Ifges(7,6) * t56;
t19 = [(Ifges(6,1) * t62 + Ifges(6,4) * t63 + t94 * t4) * t205 + t43 * t132 + (m(3) * (t89 ^ 2 + t90 ^ 2) * pkin(1) ^ 2 + t78 * mrSges(4,2) + t71 * mrSges(5,3) + Ifges(5,1) + Ifges(4,1) + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + (t223 * t95 - t225 * t92) * mrSges(6,3) + t2 * (mrSges(7,1) * t92 - t181 * t94) + t1 * (-mrSges(7,2) * t92 - t181 * t91) + t58 * (qJD(5) * t102 - t124 * t165) / 0.2e1 + t72 * (qJD(5) * t101 - t121 * t165) / 0.2e1 + t110 * t160 / 0.2e1 + t92 * t159 / 0.2e1 + t20 * t152 / 0.2e1 - 0.2e1 * mrSges(3,2) * t79 + (Ifges(6,4) * t62 + Ifges(6,2) * t63) * t237 + t62 * t129 / 0.2e1 + t71 * t134 + m(6) * (t118 * qJD(3) + (t14 * t92 + t15 * t95 + (-t34 * t92 + t35 * t95) * qJD(5)) * t70 + t117) + qJD(5) ^ 2 * t123 / 0.2e1 + m(5) * (qJD(3) * t65 + t47 * t76 + t117) + t29 * (-t106 * mrSges(7,1) - t107 * mrSges(7,2)) + m(4) * (qJD(3) * t66 - t55 * t76 + t64 * t78) + t236 * (-qJD(3) * t95 + t170 * t70) + (t161 + t43) * mrSges(5,3) + t63 * t231 + 0.2e1 * t90 * mrSges(3,1) * t176 + (t106 * t9 + t107 * t8) * mrSges(7,3) + (-m(3) * t84 - m(4) * t158 - m(5) * t142 - mrSges(2,1) * t96 - t40 * mrSges(7,1) + mrSges(2,2) * t93 - t39 * mrSges(7,2) + t228 * (-pkin(7) * t81 + t142) + t215 * t81 + (-m(7) * t139 - t213) * t82) * g(2) - t222 * t170 / 0.2e1 + (-t55 + t224) * mrSges(4,3) + (t47 + t224) * mrSges(5,2) + t221 * t67 + (t8 * mrSges(7,1) - t9 * mrSges(7,2) - t216) * t168 + m(7) * (t1 * t27 + t2 * t26 + t6 * t9 + t7 * t8) + (-m(7) * t11 - t197) * t70 * t95 + (qJD(5) * t103 - t127 * t165) * t207 + (Ifges(7,3) * t92 + t122 * t95) * t209 + (Ifges(7,6) * t92 + t125 * t95) * t210 + (Ifges(7,5) * t92 + t128 * t95) * t211 + t45 * t191 + t64 * mrSges(4,2) + qJD(4) * t60 + t11 * t51 + t6 * t32 + t7 * t33 + t26 * t12 + t27 * t13 + (m(3) * t201 + mrSges(2,1) * t93 - t38 * mrSges(7,1) + mrSges(2,2) * t96 - t37 * mrSges(7,2) + t230 * t147 + t228 * (-pkin(7) * t82 + t147) + t215 * t82 + (m(4) * pkin(2) - m(7) * t109 + t178 * t229 + t213) * t81) * g(1) + (0.2e1 * Ifges(6,5) * t205 - Ifges(6,6) * t92) * qJDD(5) + t111 * t143 + qJD(5) * t112 + t95 * t3 * t206 - (t94 * t20 + t91 * t21) * t165 / 0.2e1; (m(3) - t230) * qJDD(2) + (-m(3) - t149) * g(3) - t212 * t92 + t99 * t95; -t94 * t12 - t91 * t13 + t226 * t97 + t227 * qJDD(1) + t120 * qJD(6) + m(7) * (-t1 * t91 - t2 * t94 + (t8 * t91 - t9 * t94) * qJD(6)) + m(4) * t64 - m(6) * t43 + (-m(4) * t66 + t177 * t95 - t114 * t92 - m(7) * (t186 * t9 - t188 * t8 - t29 * t95) - m(6) * t118) * qJD(1) - t134 + (-g(1) * t81 + g(2) * t82) * t149 + (-qJD(1) * t65 - t43) * m(5); qJDD(1) * mrSges(5,2) - t97 * mrSges(5,3) + m(5) * t47 + t212 * t95 + t99 * t92 + (-m(7) * t135 - t229 * t54 + t119 - t60) * qJD(1) + (m(5) - t228) * t217; (t108 * t95 + t148 * t92 + t133) * t217 + (t108 * t92 - t148 * t95 + t132) * g(3) + (t222 / 0.2e1 + t214) * t175 + t21 * t166 / 0.2e1 + (-t110 / 0.2e1 + t111 / 0.2e1) * t97 + t11 * t131 + (t151 - t236) * t35 + (-t155 - t67) * t34 + (-pkin(5) * t11 - t17 * t8 - t18 * t9) * m(7) + (m(7) * (-qJD(6) * t135 + t136) - t33 * t166 - t32 * t167 + t218) * pkin(8) + (-t198 + t219) * mrSges(7,3) + t216 * t174 + t214 * qJD(6) + (-t112 - t8 * (mrSges(7,1) * t95 + t185 * t92) - t9 * (-mrSges(7,2) * t95 + mrSges(7,3) * t188)) * qJD(1) + t121 * t209 + t124 * t210 + t127 * t211 + t1 * t185 + t94 * t3 / 0.2e1 + t91 * t4 / 0.2e1 + Ifges(6,5) * t62 + Ifges(6,6) * t63 - t18 * t32 - t17 * t33 - t14 * mrSges(6,2) + t15 * mrSges(6,1) - pkin(5) * t5 + t123 * t143 + Ifges(6,3) * qJDD(5) + (t122 * t72 + t125 * t58 + t128 * t59) * qJD(6) / 0.2e1 - (t72 * t101 + t58 * t102 + t103 * t59) * qJD(1) / 0.2e1; -t1 * mrSges(7,2) + t2 * mrSges(7,1) - t29 * (mrSges(7,1) * t59 + mrSges(7,2) * t58) + (Ifges(7,1) * t58 - t194) * t234 + t20 * t207 + (Ifges(7,5) * t58 - Ifges(7,6) * t59) * t233 - t8 * t32 + t9 * t33 - g(1) * (mrSges(7,1) * t39 - mrSges(7,2) * t40) - g(2) * (-mrSges(7,1) * t37 + mrSges(7,2) * t38) + g(3) * t51 + (t58 * t8 + t59 * t9) * mrSges(7,3) + t159 + (-Ifges(7,2) * t59 + t21 + t52) * t235;];
tau  = t19;
