% Calculate vector of inverse dynamics joint torques for
% S5RPRPR10
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR10_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR10_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:25:55
% EndTime: 2019-12-31 18:25:59
% DurationCPUTime: 2.15s
% Computational Cost: add. (2303->247), mult. (3397->325), div. (0->0), fcn. (1654->10), ass. (0->122)
t187 = -mrSges(5,2) + mrSges(6,3);
t102 = cos(qJ(5));
t99 = sin(qJ(5));
t68 = -mrSges(6,1) * t102 + mrSges(6,2) * t99;
t186 = m(6) * pkin(4) + mrSges(5,1) - t68;
t129 = qJD(1) - qJD(3);
t169 = t99 / 0.2e1;
t152 = t99 * mrSges(6,3);
t63 = qJD(5) * mrSges(6,1) + t129 * t152;
t142 = t102 * t129;
t64 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t142;
t176 = -t102 * t64 + t99 * t63;
t185 = -t129 * mrSges(5,2) + t176;
t112 = mrSges(6,1) * t99 + mrSges(6,2) * t102;
t100 = sin(qJ(3));
t103 = cos(qJ(3));
t135 = qJ(2) * qJD(1);
t105 = -pkin(1) - pkin(2);
t77 = qJD(1) * t105 + qJD(2);
t41 = t100 * t77 + t103 * t135;
t97 = sin(pkin(8));
t162 = t41 * t97;
t40 = -t100 * t135 + t103 * t77;
t30 = -pkin(3) * t129 + t40;
t98 = cos(pkin(8));
t13 = t30 * t98 - t162;
t9 = pkin(4) * t129 - t13;
t184 = t9 * t112 + qJD(5) * (Ifges(6,5) * t102 - Ifges(6,6) * t99) / 0.2e1;
t183 = -m(4) - m(3);
t180 = mrSges(2,1) + mrSges(3,1);
t179 = -mrSges(3,3) + mrSges(2,2);
t54 = t100 * t98 + t103 * t97;
t149 = t129 * t54;
t53 = t100 * t97 - t98 * t103;
t178 = t129 * t53;
t66 = -qJ(2) * t100 + t103 * t105;
t132 = qJD(1) * qJD(2);
t79 = qJDD(1) * qJ(2) + t132;
t140 = qJD(5) * t99;
t95 = -qJDD(1) + qJDD(3);
t49 = t102 * t95 + t129 * t140;
t32 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t49;
t134 = qJD(5) * t102;
t50 = -t129 * t134 + t95 * t99;
t33 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t50;
t177 = t102 * t32 - t99 * t33;
t24 = -mrSges(6,1) * t49 + mrSges(6,2) * t50;
t73 = qJDD(1) * t105 + qJDD(2);
t21 = -qJD(3) * t41 - t100 * t79 + t103 * t73;
t12 = pkin(3) * t95 + t21;
t20 = qJD(3) * t40 + t100 * t73 + t103 * t79;
t5 = t12 * t98 - t20 * t97;
t3 = -pkin(4) * t95 - t5;
t175 = m(6) * t3 + t24;
t101 = sin(qJ(1));
t104 = cos(qJ(1));
t130 = qJ(3) + pkin(8);
t120 = sin(t130);
t121 = cos(t130);
t42 = -t101 * t120 - t104 * t121;
t43 = -t101 * t121 + t104 * t120;
t136 = t104 * t103;
t138 = t101 * t100;
t55 = -t136 - t138;
t137 = t101 * t103;
t139 = t100 * t104;
t56 = -t137 + t139;
t174 = -t55 * mrSges(4,1) - t56 * mrSges(4,2) - t186 * t42 + t187 * t43;
t173 = t56 * mrSges(4,1) - t55 * mrSges(4,2) + t186 * t43 + t187 * t42;
t44 = t68 * t129;
t172 = m(5) * t13 - m(6) * t9 - mrSges(5,1) * t129 + t44;
t34 = t98 * t41;
t14 = t97 * t30 + t34;
t10 = -pkin(7) * t129 + t14;
t7 = qJD(4) * t102 - t10 * t99;
t8 = qJD(4) * t99 + t10 * t102;
t117 = t102 * t7 + t8 * t99;
t6 = t97 * t12 + t98 * t20;
t4 = pkin(7) * t95 + t6;
t1 = qJD(5) * t7 + qJDD(4) * t99 + t102 * t4;
t163 = t1 * t102;
t2 = -qJD(5) * t8 + qJDD(4) * t102 - t4 * t99;
t108 = -qJD(5) * t117 - t2 * t99 + t163;
t171 = m(6) * t108 - t63 * t134 - t64 * t140 + t177;
t118 = t102 * t8 - t7 * t99;
t170 = m(5) * t14 + m(6) * t118 - t185;
t168 = pkin(3) * t97;
t167 = pkin(3) * t98;
t166 = pkin(7) * t43;
t164 = Ifges(6,4) * t99;
t159 = t95 * mrSges(5,1);
t158 = t95 * mrSges(5,2);
t157 = mrSges(4,1) * t129;
t155 = mrSges(4,2) * t129;
t153 = t129 * t99;
t61 = -pkin(3) + t66;
t67 = qJ(2) * t103 + t100 * t105;
t28 = t97 * t61 + t98 * t67;
t147 = t104 * pkin(1) + t101 * qJ(2);
t146 = Ifges(6,4) * t102;
t145 = Ifges(6,2) * t102;
t133 = qJDD(1) * mrSges(3,1);
t82 = pkin(3) * t139;
t125 = -pkin(7) * t42 - t82;
t93 = t104 * qJ(2);
t123 = -pkin(1) * t101 + t93;
t80 = pkin(3) * t138;
t89 = pkin(3) * t103 + pkin(2);
t122 = t104 * t89 + t147 + t80;
t119 = -pkin(3) * t136 - t80;
t27 = t61 * t98 - t67 * t97;
t111 = Ifges(6,1) * t102 - t164;
t70 = t145 + t164;
t36 = Ifges(6,6) * qJD(5) - t129 * t70;
t78 = Ifges(6,4) * t142;
t37 = -Ifges(6,1) * t153 + Ifges(6,5) * qJD(5) - t78;
t107 = t21 * mrSges(4,1) + t3 * t68 + t5 * mrSges(5,1) + t102 * (Ifges(6,4) * t50 + Ifges(6,2) * t49) / 0.2e1 - t20 * mrSges(4,2) + t49 * t70 / 0.2e1 + t50 * (Ifges(6,1) * t99 + t146) / 0.2e1 - t6 * mrSges(5,2) + (Ifges(6,1) * t50 + Ifges(6,4) * t49) * t169 + t37 * t134 / 0.2e1 - t36 * t140 / 0.2e1 - t2 * t152 + (Ifges(4,3) + Ifges(5,3)) * t95 + (-t134 * t7 - t140 * t8 + t163) * mrSges(6,3) - ((-Ifges(6,2) * t99 + t146) * t142 + t111 * t153) * qJD(5) / 0.2e1 + t184 * qJD(5) + (0.2e1 * Ifges(6,5) * t169 + Ifges(6,6) * t102) * qJDD(5);
t91 = -qJDD(1) * pkin(1) + qJDD(2);
t81 = pkin(3) * t137;
t39 = -qJD(2) * t100 - qJD(3) * t67;
t38 = qJD(2) * t103 + qJD(3) * t66;
t11 = [-t91 * mrSges(3,1) - t39 * t157 + t27 * t159 + pkin(1) * t133 + t38 * t155 - t28 * t158 + m(4) * (t20 * t67 + t21 * t66 + t38 * t41 + t39 * t40) + m(5) * (t27 * t5 + t28 * t6) - t107 + m(3) * (-pkin(1) * t91 + (t79 + t132) * qJ(2)) + (t66 * mrSges(4,1) - t67 * mrSges(4,2)) * t95 + t175 * (pkin(4) - t27) - t172 * (t38 * t97 - t98 * t39) + (Ifges(3,2) + Ifges(2,3)) * qJDD(1) + t170 * (t38 * t98 + t39 * t97) + 0.2e1 * t79 * mrSges(3,3) + t171 * (-pkin(7) + t28) + (-m(6) * (t122 + t166) - m(5) * t122 + t183 * t147 + (-m(4) * pkin(2) - t180) * t104 + t179 * t101 - t174) * g(2) + (-m(5) * (t82 + t93) - m(4) * t93 - m(3) * t123 - m(6) * (t123 - t125) + t179 * t104 + (-m(5) * (-pkin(1) - t89) - m(4) * t105 + m(6) * t89 + t180) * t101 - t173) * g(1); -t133 + t53 * t24 + t149 * t44 + (-m(3) * qJ(2) - mrSges(3,3)) * qJD(1) ^ 2 + (mrSges(4,1) * t103 - mrSges(5,1) * t53 - mrSges(4,2) * t100) * t95 + (-t158 + (-t102 * t63 - t99 * t64) * qJD(5) + t177) * t54 - (mrSges(5,1) * t149 + (mrSges(4,1) * t100 + mrSges(4,2) * t103) * t129) * t129 + m(3) * t91 - t185 * t178 + (-g(1) * t101 + g(2) * t104) * (m(5) + m(6) - t183) + (t108 * t54 + t118 * t178 - t149 * t9 + t3 * t53) * m(6) + (t13 * t149 + t14 * t178 - t5 * t53 + t54 * t6) * m(5) + (t100 * t20 + t103 * t21 - t129 * (-t100 * t40 + t103 * t41)) * m(4); m(5) * (t5 * t98 + t6 * t97) * pkin(3) - t158 * t168 - t40 * t155 - t41 * t157 + t159 * t167 + t107 + t175 * (-pkin(4) - t167) + t172 * (t40 * t97 + t34) - t170 * (t40 * t98 - t162) + (-m(5) * t119 - m(6) * (t119 - t166) + t174) * g(2) + (-m(5) * (-t82 + t81) - m(6) * (t125 + t81) + t173) * g(1) + t171 * (pkin(7) + t168); t102 * t33 + t99 * t32 - t176 * qJD(5) + (qJD(5) * t118 + t1 * t99 + t102 * t2 + g(3)) * m(6) + (qJDD(4) + g(3)) * m(5); t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t50 + Ifges(6,6) * t49 + Ifges(6,3) * qJDD(5) - g(3) * t68 + t8 * t63 - t7 * t64 - (t36 * t169 - (-t99 * t111 / 0.2e1 + t145 * t169) * t129 + t117 * mrSges(6,3) - (t37 - t78) * t102 / 0.2e1 - t184) * t129 + (-g(1) * t42 - g(2) * t43) * t112;];
tau = t11;
