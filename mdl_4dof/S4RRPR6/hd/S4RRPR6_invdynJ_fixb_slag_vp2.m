% Calculate vector of inverse dynamics joint torques for
% S4RRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:05
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:04:23
% EndTime: 2019-12-31 17:04:32
% DurationCPUTime: 4.70s
% Computational Cost: add. (2250->333), mult. (5399->466), div. (0->0), fcn. (3712->12), ass. (0->160)
t131 = qJDD(2) + qJDD(4);
t133 = qJD(2) + qJD(4);
t138 = sin(qJ(4));
t141 = cos(qJ(4));
t135 = sin(pkin(7));
t136 = cos(pkin(7));
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t100 = -t135 * t139 + t136 * t142;
t87 = t100 * qJD(1);
t171 = t142 * qJD(1);
t172 = t139 * qJD(1);
t88 = -t135 * t171 - t136 * t172;
t159 = t138 * t88 + t141 * t87;
t170 = qJD(1) * qJD(2);
t161 = t139 * t170;
t169 = qJDD(1) * t142;
t108 = -t161 + t169;
t109 = qJDD(1) * t139 + t142 * t170;
t63 = t108 * t136 - t109 * t135;
t64 = t108 * t135 + t109 * t136;
t14 = qJD(4) * t159 + t138 * t63 + t141 * t64;
t49 = t138 * t87 - t141 * t88;
t15 = -qJD(4) * t49 - t138 * t64 + t141 * t63;
t193 = Ifges(5,4) * t49;
t99 = t109 * pkin(5);
t53 = qJDD(2) * pkin(2) - qJ(3) * t109 - qJD(3) * t172 - t99;
t125 = pkin(5) * t169;
t175 = qJD(2) * t139;
t167 = pkin(5) * t175;
t173 = qJD(3) * t142;
t59 = qJ(3) * t108 + t125 + (-t167 + t173) * qJD(1);
t24 = -t135 * t59 + t136 * t53;
t10 = qJDD(2) * pkin(3) - pkin(6) * t64 + t24;
t25 = t135 * t53 + t136 * t59;
t16 = pkin(6) * t63 + t25;
t198 = pkin(6) * t88;
t137 = -qJ(3) - pkin(5);
t116 = t137 * t142;
t107 = qJD(1) * t116;
t91 = t135 * t107;
t114 = t137 * t139;
t106 = qJD(1) * t114;
t97 = qJD(2) * pkin(2) + t106;
t55 = t136 * t97 + t91;
t30 = qJD(2) * pkin(3) + t198 + t55;
t199 = pkin(6) * t87;
t176 = t136 * t107;
t56 = t135 * t97 - t176;
t31 = t56 + t199;
t8 = -t138 * t31 + t141 * t30;
t2 = qJD(4) * t8 + t10 * t138 + t141 * t16;
t41 = Ifges(5,4) * t159;
t22 = Ifges(5,1) * t49 + Ifges(5,5) * t133 + t41;
t9 = t138 * t30 + t141 * t31;
t3 = -qJD(4) * t9 + t10 * t141 - t138 * t16;
t190 = pkin(2) * t142;
t124 = pkin(1) + t190;
t111 = -qJD(1) * t124 + qJD(3);
t65 = -pkin(3) * t87 + t111;
t220 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t14 + Ifges(5,6) * t15 + Ifges(5,3) * t131 - (Ifges(5,5) * t159 - Ifges(5,6) * t49) * t133 / 0.2e1 + (t159 * t8 + t49 * t9) * mrSges(5,3) - (-Ifges(5,2) * t49 + t22 + t41) * t159 / 0.2e1 - t65 * (mrSges(5,1) * t49 + mrSges(5,2) * t159) - (Ifges(5,1) * t159 - t193) * t49 / 0.2e1;
t61 = -t106 * t135 + t176;
t33 = t61 - t199;
t62 = t136 * t106 + t91;
t34 = t62 + t198;
t123 = pkin(2) * t136 + pkin(3);
t192 = pkin(2) * t135;
t83 = t123 * t141 - t138 * t192;
t219 = t83 * qJD(4) - t138 * t33 - t141 * t34;
t84 = t123 * t138 + t141 * t192;
t218 = -t84 * qJD(4) + t138 * t34 - t141 * t33;
t140 = sin(qJ(1));
t143 = cos(qJ(1));
t217 = g(1) * t143 + g(2) * t140;
t115 = -t142 * mrSges(3,1) + t139 * mrSges(3,2);
t134 = qJ(2) + pkin(7);
t128 = sin(t134);
t129 = cos(t134);
t130 = qJ(4) + t134;
t121 = sin(t130);
t122 = cos(t130);
t156 = t122 * mrSges(5,1) - t121 * mrSges(5,2);
t216 = t129 * mrSges(4,1) - t128 * mrSges(4,2) - t115 + t156;
t21 = Ifges(5,2) * t159 + Ifges(5,6) * t133 + t193;
t214 = t21 / 0.2e1;
t178 = qJDD(2) / 0.2e1;
t177 = qJDD(1) * pkin(1);
t98 = -pkin(5) * t161 + t125;
t209 = t139 * t99 + t142 * t98;
t208 = 0.2e1 * t178;
t160 = pkin(3) * t129 + t190;
t207 = mrSges(2,1) + m(5) * (pkin(1) + t160) + m(4) * t124 + m(3) * pkin(1) + t216;
t206 = mrSges(2,2) + m(5) * (-pkin(6) + t137) - mrSges(5,3) + m(4) * t137 - mrSges(4,3) - m(3) * pkin(5) - mrSges(3,3);
t205 = m(4) * pkin(2);
t202 = t49 / 0.2e1;
t200 = -t88 / 0.2e1;
t196 = t139 / 0.2e1;
t194 = Ifges(4,4) * t88;
t191 = pkin(2) * t139;
t187 = -qJD(2) / 0.2e1;
t158 = qJD(2) * t137;
t85 = t139 * t158 + t173;
t86 = -qJD(3) * t139 + t142 * t158;
t38 = t135 * t86 + t136 * t85;
t186 = mrSges(3,2) * t142;
t185 = Ifges(3,4) * t139;
t184 = Ifges(3,4) * t142;
t183 = Ifges(3,2) * t142;
t67 = t135 * t114 - t136 * t116;
t174 = qJD(2) * t142;
t168 = pkin(2) * t175;
t163 = -t63 * mrSges(4,1) + t64 * mrSges(4,2);
t162 = -t15 * mrSges(5,1) + t14 * mrSges(5,2);
t37 = -t135 * t85 + t136 * t86;
t66 = t136 * t114 + t116 * t135;
t155 = -g(1) * t140 + g(2) * t143;
t153 = mrSges(5,1) * t121 + mrSges(5,2) * t122;
t152 = t183 + t185;
t151 = Ifges(3,5) * t142 - Ifges(3,6) * t139;
t101 = t135 * t142 + t136 * t139;
t39 = -pkin(6) * t101 + t66;
t40 = pkin(6) * t100 + t67;
t19 = -t138 * t40 + t141 * t39;
t20 = t138 * t39 + t141 * t40;
t57 = t100 * t141 - t101 * t138;
t58 = t100 * t138 + t101 * t141;
t149 = pkin(1) * (mrSges(3,1) * t139 + t186);
t81 = -pkin(2) * t108 + qJDD(3) - t177;
t148 = t139 * (Ifges(3,1) * t142 - t185);
t126 = Ifges(3,4) * t171;
t113 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t171;
t112 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t172;
t96 = Ifges(3,1) * t172 + Ifges(3,5) * qJD(2) + t126;
t95 = Ifges(3,6) * qJD(2) + qJD(1) * t152;
t90 = t100 * qJD(2);
t89 = t101 * qJD(2);
t82 = Ifges(4,4) * t87;
t72 = -pkin(3) * t100 - t124;
t71 = qJD(2) * mrSges(4,1) + mrSges(4,3) * t88;
t70 = -qJD(2) * mrSges(4,2) + mrSges(4,3) * t87;
t69 = pkin(3) * t89 + t168;
t68 = pkin(2) * t172 - pkin(3) * t88;
t54 = -mrSges(4,1) * t87 - mrSges(4,2) * t88;
t51 = qJDD(2) * mrSges(4,1) - mrSges(4,3) * t64;
t50 = -qJDD(2) * mrSges(4,2) + mrSges(4,3) * t63;
t45 = -t88 * Ifges(4,1) + Ifges(4,5) * qJD(2) + t82;
t44 = t87 * Ifges(4,2) + Ifges(4,6) * qJD(2) - t194;
t36 = mrSges(5,1) * t133 - mrSges(5,3) * t49;
t35 = -mrSges(5,2) * t133 + mrSges(5,3) * t159;
t32 = -pkin(3) * t63 + t81;
t29 = -pkin(6) * t89 + t38;
t28 = -pkin(6) * t90 + t37;
t27 = -qJD(4) * t58 - t138 * t90 - t141 * t89;
t26 = qJD(4) * t57 - t138 * t89 + t141 * t90;
t23 = -mrSges(5,1) * t159 + mrSges(5,2) * t49;
t7 = -mrSges(5,2) * t131 + mrSges(5,3) * t15;
t6 = mrSges(5,1) * t131 - mrSges(5,3) * t14;
t5 = -qJD(4) * t20 - t138 * t29 + t141 * t28;
t4 = qJD(4) * t19 + t138 * t28 + t141 * t29;
t1 = [qJD(2) ^ 2 * t151 / 0.2e1 + t108 * t152 / 0.2e1 + (Ifges(3,1) * t109 + Ifges(3,4) * t108 + Ifges(3,5) * qJDD(2)) * t196 + (Ifges(5,1) * t26 + Ifges(5,4) * t27) * t202 + (Ifges(3,5) * t139 + Ifges(3,6) * t142) * t178 + t109 * (t139 * Ifges(3,1) + t184) / 0.2e1 + (mrSges(5,2) * t32 - mrSges(5,3) * t3 + Ifges(5,1) * t14 + Ifges(5,4) * t15 + Ifges(5,5) * t131) * t58 + (-t26 * t8 + t27 * t9) * mrSges(5,3) + t142 * (Ifges(3,4) * t109 + Ifges(3,2) * t108 + Ifges(3,6) * qJDD(2)) / 0.2e1 + (t148 + t142 * (-Ifges(3,2) * t139 + t184)) * t170 / 0.2e1 + t159 * (Ifges(5,4) * t26 + Ifges(5,2) * t27) / 0.2e1 + t133 * (Ifges(5,5) * t26 + Ifges(5,6) * t27) / 0.2e1 - t89 * t44 / 0.2e1 + t90 * t45 / 0.2e1 + t38 * t70 + t37 * t71 + t65 * (-mrSges(5,1) * t27 + mrSges(5,2) * t26) + t66 * t51 + t67 * t50 + t69 * t23 + t4 * t35 + t5 * t36 + t26 * t22 / 0.2e1 + t19 * t6 + t20 * t7 + t27 * t214 + t54 * t168 + (Ifges(4,1) * t90 - Ifges(4,4) * t89) * t200 + (-t55 * t90 - t56 * t89) * mrSges(4,3) + t111 * (mrSges(4,1) * t89 + mrSges(4,2) * t90) + t87 * (Ifges(4,4) * t90 - Ifges(4,2) * t89) / 0.2e1 + qJD(2) * (Ifges(4,5) * t90 - Ifges(4,6) * t89) / 0.2e1 + (m(3) * t177 + mrSges(3,1) * t108 - mrSges(3,2) * t109) * pkin(1) + (-t81 * mrSges(4,1) + t25 * mrSges(4,3) + Ifges(4,4) * t64 + t63 * Ifges(4,2) + t208 * Ifges(4,6)) * t100 + (t81 * mrSges(4,2) - t24 * mrSges(4,3) + t64 * Ifges(4,1) + Ifges(4,4) * t63 + Ifges(4,5) * t208) * t101 + t209 * mrSges(3,3) + (t142 * (-qJDD(2) * mrSges(3,2) + mrSges(3,3) * t108) - t139 * (qJDD(2) * mrSges(3,1) - mrSges(3,3) * t109) - t112 * t174 + m(3) * t209) * pkin(5) + (t140 * t207 + t143 * t206) * g(1) + (t140 * t206 - t143 * t207) * g(2) + m(5) * (t19 * t3 + t2 * t20 + t32 * t72 + t4 * t9 + t5 * t8 + t65 * t69) - t113 * t167 + m(4) * (t111 * t168 - t124 * t81 + t24 * t66 + t25 * t67 + t37 * t55 + t38 * t56) - t149 * t170 - t115 * t177 + t96 * t174 / 0.2e1 - t95 * t175 / 0.2e1 + t72 * t162 - t124 * t163 + (-mrSges(5,1) * t32 + mrSges(5,3) * t2 + Ifges(5,4) * t14 + Ifges(5,2) * t15 + Ifges(5,6) * t131) * t57 + Ifges(2,3) * qJDD(1); (Ifges(4,3) + Ifges(3,3)) * qJDD(2) + (t135 * t50 + t136 * t51) * pkin(2) - m(4) * (t55 * t61 + t56 * t62) + (t55 * t87 - t56 * t88) * mrSges(4,3) + t88 * (Ifges(4,1) * t87 + t194) / 0.2e1 + t44 * t200 + (t135 * t25 + t136 * t24) * t205 + (Ifges(4,5) * t87 + Ifges(4,6) * t88) * t187 + t218 * t36 + (-g(3) * t160 + t2 * t84 + t218 * t8 + t219 * t9 + t3 * t83 - t65 * t68) * m(5) + t219 * t35 - (Ifges(4,2) * t88 + t45 + t82) * t87 / 0.2e1 + t49 * t214 + (-m(4) * t190 - t216) * g(3) + (t151 * t187 + t95 * t196 + (t183 * t196 - t148 / 0.2e1 + t149) * qJD(1) + (t142 * t112 + t139 * t113) * pkin(5) + (-m(4) * t111 - t54) * t191 - (t126 + t96) * t142 / 0.2e1) * qJD(1) + Ifges(3,6) * t108 + Ifges(3,5) * t109 - t111 * (-mrSges(4,1) * t88 + mrSges(4,2) * t87) - t98 * mrSges(3,2) - t99 * mrSges(3,1) + t83 * t6 + t84 * t7 - t62 * t70 - t61 * t71 + Ifges(4,6) * t63 + Ifges(4,5) * t64 - t68 * t23 - t25 * mrSges(4,2) + t24 * mrSges(4,1) + t217 * (-m(5) * (-pkin(3) * t128 - t191) + mrSges(4,1) * t128 + t186 + mrSges(4,2) * t129 + (mrSges(3,1) + t205) * t139 + t153) + t220; -t159 * t35 + t49 * t36 - t87 * t70 - t88 * t71 + t162 + t163 + (-t159 * t9 + t49 * t8 + t155 + t32) * m(5) + (-t55 * t88 - t56 * t87 + t155 + t81) * m(4); -g(3) * t156 + t217 * t153 + t21 * t202 - t8 * t35 + t9 * t36 + t220;];
tau = t1;
