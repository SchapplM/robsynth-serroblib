% Calculate vector of inverse dynamics joint torques for
% S4PRRR6
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
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
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
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:34:40
% EndTime: 2019-12-31 16:34:45
% DurationCPUTime: 2.78s
% Computational Cost: add. (1168->256), mult. (2678->377), div. (0->0), fcn. (1689->10), ass. (0->126)
t94 = sin(qJ(3));
t97 = cos(qJ(3));
t120 = -mrSges(4,1) * t97 + mrSges(4,2) * t94;
t80 = pkin(3) * t97 + pkin(2);
t90 = qJ(3) + qJ(4);
t84 = sin(t90);
t85 = cos(t90);
t176 = m(4) * pkin(2) + m(5) * t80 + t85 * mrSges(5,1) - t84 * mrSges(5,2) + mrSges(3,1) - t120;
t99 = -pkin(6) - pkin(5);
t175 = -m(4) * pkin(5) + m(5) * t99 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t130 = qJD(2) * qJD(3);
t68 = qJDD(2) * t97 - t130 * t94;
t174 = t68 / 0.2e1;
t93 = sin(qJ(4));
t96 = cos(qJ(4));
t114 = t93 * t94 - t96 * t97;
t98 = cos(qJ(2));
t136 = qJD(1) * t98;
t74 = t99 * t94;
t75 = t99 * t97;
t37 = t74 * t96 + t75 * t93;
t125 = qJD(3) * t99;
t66 = t94 * t125;
t67 = t97 * t125;
t173 = qJD(4) * t37 + t114 * t136 + t66 * t96 + t67 * t93;
t64 = t93 * t97 + t94 * t96;
t109 = t64 * t98;
t38 = t74 * t93 - t75 * t96;
t172 = qJD(1) * t109 - qJD(4) * t38 - t66 * t93 + t67 * t96;
t58 = t114 * qJD(2);
t59 = t64 * qJD(2);
t30 = mrSges(5,1) * t58 + mrSges(5,2) * t59;
t171 = -t120 * qJD(2) - t30;
t95 = sin(qJ(2));
t52 = t114 * t95;
t165 = m(5) * pkin(3);
t170 = -t165 - mrSges(4,1);
t133 = t94 * qJD(2);
t72 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t133;
t132 = t97 * qJD(2);
t73 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t132;
t169 = t97 * t72 + t94 * t73;
t69 = qJDD(2) * t94 + t130 * t97;
t168 = t97 * (-qJDD(3) * mrSges(4,2) + mrSges(4,3) * t68) - t94 * (qJDD(3) * mrSges(4,1) - mrSges(4,3) * t69);
t167 = t72 * t94 - t73 * t97;
t135 = qJD(3) * t94;
t131 = qJD(1) * qJD(2);
t79 = t98 * t131;
t71 = t95 * qJDD(1) + t79;
t61 = qJDD(2) * pkin(5) + t71;
t137 = qJD(1) * t95;
t76 = qJD(2) * pkin(5) + t137;
t31 = -t135 * t76 + t97 * t61;
t134 = qJD(3) * t97;
t32 = -t134 * t76 - t61 * t94;
t115 = t31 * t97 - t32 * t94;
t87 = qJD(3) + qJD(4);
t100 = qJD(2) ^ 2;
t163 = t59 / 0.2e1;
t160 = g(3) * t95;
t159 = Ifges(4,4) * t94;
t158 = Ifges(4,4) * t97;
t123 = pkin(6) * qJD(2) + t76;
t50 = t123 * t97;
t154 = t50 * t93;
t49 = t123 * t94;
t45 = qJD(3) * pkin(3) - t49;
t21 = t45 * t96 - t154;
t157 = t21 * mrSges(5,3);
t153 = t50 * t96;
t152 = t59 * mrSges(5,3);
t151 = t59 * Ifges(5,4);
t91 = sin(pkin(7));
t148 = t91 * t98;
t92 = cos(pkin(7));
t147 = t92 * t98;
t144 = t94 * t98;
t141 = t97 * t98;
t140 = (-t148 * t84 - t85 * t92) * mrSges(5,1) + (-t148 * t85 + t84 * t92) * mrSges(5,2);
t139 = (-t147 * t84 + t85 * t91) * mrSges(5,1) + (-t147 * t85 - t84 * t91) * mrSges(5,2);
t129 = pkin(3) * t133;
t128 = pkin(3) * t135;
t78 = t95 * t131;
t70 = qJDD(1) * t98 - t78;
t119 = mrSges(4,1) * t94 + mrSges(4,2) * t97;
t118 = -mrSges(5,1) * t84 - mrSges(5,2) * t85;
t117 = Ifges(4,2) * t97 + t159;
t116 = Ifges(4,5) * t97 - Ifges(4,6) * t94;
t22 = t45 * t93 + t153;
t77 = -qJD(2) * pkin(2) - t136;
t111 = t77 * t119;
t110 = t94 * (Ifges(4,1) * t97 - t159);
t60 = -qJDD(2) * pkin(2) - t70;
t108 = t64 * qJD(4);
t107 = t114 * qJD(4);
t104 = t77 * t95 + (t94 ^ 2 + t97 ^ 2) * t98 * t76;
t19 = -qJD(2) * t107 + t68 * t93 + t69 * t96;
t18 = qJDD(3) * pkin(3) - pkin(6) * t69 + t32;
t23 = pkin(6) * t68 + t31;
t2 = qJD(4) * t21 + t18 * t93 + t23 * t96;
t20 = -qJD(2) * t108 + t68 * t96 - t69 * t93;
t27 = -t58 * Ifges(5,2) + t87 * Ifges(5,6) + t151;
t55 = Ifges(5,4) * t58;
t28 = t59 * Ifges(5,1) + t87 * Ifges(5,5) - t55;
t3 = -qJD(4) * t22 + t18 * t96 - t23 * t93;
t62 = -qJD(2) * t80 - t136;
t86 = qJDD(3) + qJDD(4);
t103 = t3 * mrSges(5,1) - t2 * mrSges(5,2) - t58 * t157 + t27 * t163 - t62 * (mrSges(5,1) * t59 - mrSges(5,2) * t58) - t59 * (-Ifges(5,1) * t58 - t151) / 0.2e1 + Ifges(5,6) * t20 + Ifges(5,5) * t19 - t87 * (-Ifges(5,5) * t58 - Ifges(5,6) * t59) / 0.2e1 + Ifges(5,3) * t86 + (-Ifges(5,2) * t59 + t28 - t55) * t58 / 0.2e1;
t34 = -qJD(3) * t64 - t108;
t81 = Ifges(4,4) * t132;
t57 = Ifges(4,1) * t133 + Ifges(4,5) * qJD(3) + t81;
t56 = Ifges(4,6) * qJD(3) + qJD(2) * t117;
t51 = t64 * t95;
t40 = mrSges(5,1) * t87 - t152;
t39 = -mrSges(5,2) * t87 - mrSges(5,3) * t58;
t36 = -mrSges(4,1) * t68 + mrSges(4,2) * t69;
t35 = -pkin(3) * t68 + t60;
t33 = -qJD(3) * t114 - t107;
t25 = -t49 * t96 - t154;
t24 = t49 * t93 - t153;
t13 = -mrSges(5,2) * t86 + mrSges(5,3) * t20;
t12 = mrSges(5,1) * t86 - mrSges(5,3) * t19;
t7 = -qJD(2) * t109 + t52 * t87;
t6 = t34 * t95 - t58 * t98;
t4 = -mrSges(5,1) * t20 + mrSges(5,2) * t19;
t1 = [m(2) * qJDD(1) - t51 * t12 - t52 * t13 + t6 * t39 + t7 * t40 + (qJDD(2) * mrSges(3,1) - t100 * mrSges(3,2) - qJD(2) * t167 - t36 - t4) * t98 + (-t100 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - qJD(2) * t171 - qJD(3) * t169 + t168) * t95 + m(4) * (qJD(2) * t104 + t115 * t95 - t60 * t98) + m(5) * (qJD(2) * t62 * t95 - t2 * t52 + t21 * t7 + t22 * t6 - t3 * t51 - t35 * t98) + m(3) * (t70 * t98 + t71 * t95) + (-m(2) - m(3) - m(4) - m(5)) * g(3); t60 * t120 + (-pkin(2) * t60 - t104 * qJD(1)) * m(4) + t172 * t40 + (t2 * t38 + t3 * t37 - t35 * t80 + (t128 - t137) * t62 + t173 * t22 + t172 * t21) * m(5) + t173 * t39 + (Ifges(4,1) * t69 + Ifges(4,4) * t174) * t94 + t97 * (Ifges(4,4) * t69 + Ifges(4,2) * t68) / 0.2e1 - (-mrSges(5,1) * t35 + mrSges(5,3) * t2 + Ifges(5,4) * t19 + Ifges(5,2) * t20 + Ifges(5,6) * t86) * t114 + (t97 * (-Ifges(4,2) * t94 + t158) + t110) * t130 / 0.2e1 + (mrSges(5,2) * t35 - mrSges(5,3) * t3 + Ifges(5,1) * t19 + Ifges(5,4) * t20 + Ifges(5,5) * t86) * t64 + (t111 + t116 * qJD(3) / 0.2e1) * qJD(3) + (-t71 + t79) * mrSges(3,2) + (t70 + t78) * mrSges(3,1) + (Ifges(5,1) * t33 + Ifges(5,4) * t34) * t163 + t117 * t174 + t171 * t137 + (m(4) * t115 - t134 * t72 - t135 * t73 + t168) * pkin(5) + qJDD(3) * (Ifges(4,5) * t94 + Ifges(4,6) * t97) + t87 * (Ifges(5,5) * t33 + Ifges(5,6) * t34) / 0.2e1 - t80 * t4 - t58 * (Ifges(5,4) * t33 + Ifges(5,2) * t34) / 0.2e1 + t62 * (-mrSges(5,1) * t34 + mrSges(5,2) * t33) + t37 * t12 + t38 * t13 + t33 * t28 / 0.2e1 + t34 * t27 / 0.2e1 - pkin(2) * t36 + t22 * t34 * mrSges(5,3) + (g(1) * t92 + g(2) * t91) * (t175 * t98 + t176 * t95) + (t175 * t95 - t176 * t98) * g(3) + t69 * t158 / 0.2e1 + t30 * t128 - t33 * t157 + t57 * t134 / 0.2e1 - t56 * t135 / 0.2e1 + t115 * mrSges(4,3) + t167 * t136 + Ifges(3,3) * qJDD(2); -qJD(2) * t111 + (t2 * t93 + t3 * t96 + (-t21 * t93 + t22 * t96) * qJD(4)) * t165 + t22 * t152 + Ifges(4,5) * t69 + Ifges(4,6) * t68 - t25 * t39 - t24 * t40 - t31 * mrSges(4,2) + t32 * mrSges(4,1) + t103 - t100 * t110 / 0.2e1 + t56 * t133 / 0.2e1 - t116 * t130 / 0.2e1 - t30 * t129 - m(5) * (t129 * t62 + t21 * t24 + t22 * t25) + Ifges(4,3) * qJDD(3) + t169 * t76 + (t165 * t94 - t118 + t119) * t160 - (-Ifges(4,2) * t133 + t57 + t81) * t132 / 0.2e1 + (-t140 - (-t141 * t91 + t92 * t94) * mrSges(4,2) + t170 * (-t144 * t91 - t92 * t97)) * g(2) + (-t139 - (-t141 * t92 - t91 * t94) * mrSges(4,2) + t170 * (-t144 * t92 + t91 * t97)) * g(1) + ((t39 * t96 - t40 * t93) * qJD(4) + t96 * t12 + t93 * t13) * pkin(3); (t40 + t152) * t22 - t118 * t160 - t21 * t39 - g(2) * t140 - g(1) * t139 + t103;];
tau = t1;
