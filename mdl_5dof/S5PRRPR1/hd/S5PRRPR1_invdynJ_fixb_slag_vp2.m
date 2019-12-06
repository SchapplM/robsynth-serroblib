% Calculate vector of inverse dynamics joint torques for
% S5PRRPR1
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:16
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:15:34
% EndTime: 2019-12-05 16:15:39
% DurationCPUTime: 1.44s
% Computational Cost: add. (1750->225), mult. (2640->303), div. (0->0), fcn. (1601->12), ass. (0->111)
t127 = cos(qJ(3));
t160 = pkin(2) * qJD(2);
t137 = -t127 * t160 + qJD(4);
t124 = sin(qJ(5));
t126 = cos(qJ(5));
t121 = sin(pkin(9));
t123 = -pkin(7) - qJ(4);
t82 = t123 * t121;
t122 = cos(pkin(9));
t113 = t122 * pkin(7);
t83 = qJ(4) * t122 + t113;
t44 = -t124 * t83 + t126 * t82;
t72 = -t121 * t124 + t122 * t126;
t176 = t44 * qJD(5) + t137 * t72;
t45 = t124 * t82 + t126 * t83;
t73 = t121 * t126 + t122 * t124;
t175 = -t45 * qJD(5) - t137 * t73;
t136 = -t122 * mrSges(5,1) + t121 * mrSges(5,2);
t174 = mrSges(4,1) - t136;
t125 = sin(qJ(3));
t173 = t125 * t160;
t118 = pkin(9) + qJ(5);
t108 = sin(t118);
t110 = cos(t118);
t172 = t110 * mrSges(6,1) - t108 * mrSges(6,2);
t171 = mrSges(4,2) - mrSges(6,3) - mrSges(5,3);
t170 = -t172 - t174;
t169 = -m(2) - m(3);
t120 = qJD(2) + qJD(3);
t147 = t121 ^ 2 + t122 ^ 2;
t168 = t120 * t147;
t107 = t122 * qJD(1);
t78 = qJ(4) * t120 + t173;
t54 = -t121 * t78 + t107;
t55 = t121 * qJD(1) + t122 * t78;
t134 = -t121 * t54 + t122 * t55;
t58 = t72 * t120;
t59 = t73 * t120;
t167 = (mrSges(6,1) * t58 - mrSges(6,2) * t59 + t174 * t120) * t125 + mrSges(4,2) * t120 * t127;
t166 = m(4) * pkin(2);
t164 = t59 / 0.2e1;
t163 = Ifges(6,4) * t59;
t162 = pkin(2) * t127;
t119 = pkin(8) + qJ(2);
t112 = qJ(3) + t119;
t97 = sin(t112);
t98 = cos(t112);
t161 = t98 * pkin(3) + t97 * qJ(4);
t159 = pkin(2) * qJD(3);
t105 = t122 * qJDD(1);
t115 = qJDD(2) + qJDD(3);
t142 = qJD(2) * t159;
t150 = pkin(2) * qJDD(2);
t71 = t125 * t150 + t127 * t142;
t48 = qJ(4) * t115 + qJD(4) * t120 + t71;
t36 = -t121 * t48 + t105;
t155 = t121 * t36;
t37 = t121 * qJDD(1) + t122 * t48;
t152 = t122 * t37;
t149 = t115 * t121;
t148 = t115 * t122;
t146 = m(4) + m(5) + m(6);
t143 = t125 * t159;
t100 = t122 * pkin(4) + pkin(3);
t64 = t72 * qJD(5);
t28 = t73 * t115 + t120 * t64;
t65 = t73 * qJD(5);
t29 = t72 * t115 - t120 * t65;
t6 = -t29 * mrSges(6,1) + t28 * mrSges(6,2);
t141 = t98 * t100 - t123 * t97;
t140 = t147 * t115;
t62 = -mrSges(5,1) * t148 + mrSges(5,2) * t149;
t138 = -g(1) * t97 + g(2) * t98;
t70 = -t125 * t142 + t127 * t150;
t40 = t107 + (-pkin(7) * t120 - t78) * t121;
t41 = t120 * t113 + t55;
t14 = -t124 * t41 + t126 * t40;
t15 = t124 * t40 + t126 * t41;
t99 = pkin(2) * t125 + qJ(4);
t68 = (-pkin(7) - t99) * t121;
t69 = t122 * t99 + t113;
t34 = -t124 * t69 + t126 * t68;
t35 = t124 * t68 + t126 * t69;
t133 = qJDD(4) - t70;
t130 = t170 * t98 + t171 * t97;
t129 = (m(5) * pkin(3) + m(6) * t100 - t170) * t97 + (-m(5) * qJ(4) + m(6) * t123 + t171) * t98;
t30 = t105 + (-pkin(7) * t115 - t48) * t121;
t31 = pkin(7) * t148 + t37;
t2 = t14 * qJD(5) + t124 * t30 + t126 * t31;
t22 = Ifges(6,2) * t58 + Ifges(6,6) * qJD(5) + t163;
t53 = Ifges(6,4) * t58;
t23 = Ifges(6,1) * t59 + Ifges(6,5) * qJD(5) + t53;
t3 = -t15 * qJD(5) - t124 * t31 + t126 * t30;
t42 = -t100 * t115 + t133;
t56 = -t100 * t120 + t137;
t57 = -t115 * pkin(3) + t133;
t128 = Ifges(4,3) * t115 + t64 * t23 / 0.2e1 - t65 * t22 / 0.2e1 + t58 * (Ifges(6,4) * t64 - Ifges(6,2) * t65) / 0.2e1 + (Ifges(6,1) * t64 - Ifges(6,4) * t65) * t164 + t56 * (mrSges(6,1) * t65 + mrSges(6,2) * t64) + qJD(5) * (Ifges(6,5) * t64 - Ifges(6,6) * t65) / 0.2e1 + mrSges(5,3) * t152 + t57 * t136 + (Ifges(5,1) * t121 + Ifges(5,4) * t122) * t149 + (Ifges(5,4) * t121 + Ifges(5,2) * t122) * t148 + t70 * mrSges(4,1) - t71 * mrSges(4,2) + (-t14 * t64 - t15 * t65) * mrSges(6,3) + (t42 * mrSges(6,2) - t3 * mrSges(6,3) + Ifges(6,1) * t28 + Ifges(6,4) * t29 + Ifges(6,5) * qJDD(5)) * t73 + (-t42 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t28 + Ifges(6,2) * t29 + Ifges(6,6) * qJDD(5)) * t72;
t111 = cos(t119);
t109 = sin(t119);
t101 = -pkin(3) - t162;
t96 = pkin(2) * t111;
t92 = t127 * t159 + qJD(4);
t81 = -t100 - t162;
t75 = -t120 * pkin(3) + t137;
t47 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t59;
t46 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t58;
t21 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t29;
t20 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t28;
t11 = -t35 * qJD(5) - t73 * t92;
t10 = t34 * qJD(5) + t72 * t92;
t1 = [t72 * t20 + t73 * t21 + t64 * t46 - t65 * t47 + m(5) * (t121 * t37 + t122 * t36) + m(6) * (-t14 * t65 + t15 * t64 + t2 * t73 + t3 * t72) + (m(4) - t169) * qJDD(1) + (-t146 + t169) * g(3); m(6) * (t10 * t15 + t11 * t14 + t56 * t143 + t2 * t35 + t3 * t34 + t42 * t81) + (t99 * t140 + t92 * t168 - t155) * mrSges(5,3) + (t111 * mrSges(3,2) + (t146 * pkin(2) + mrSges(3,1)) * t109 + t129) * g(1) + (-m(6) * (t141 + t96) - m(5) * (t96 + t161) + t109 * mrSges(3,2) + (-mrSges(3,1) - t166) * t111 + t130) * g(2) + ((mrSges(4,1) * t127 - mrSges(4,2) * t125) * t115 - t167 * qJD(3)) * pkin(2) + m(5) * (t75 * t143 + t101 * t57 + (t37 * t99 + t55 * t92) * t122 + (-t36 * t99 - t54 * t92) * t121) + t101 * t62 + t81 * t6 + t10 * t46 + t11 * t47 + t34 * t20 + t35 * t21 + t128 + (t125 * t71 + t127 * t70) * t166 + Ifges(3,3) * qJDD(2); (qJ(4) * t140 + t137 * t168 - t155) * mrSges(5,3) + t167 * t160 + t129 * g(1) + t130 * g(2) + t175 * t47 + t176 * t46 - t100 * t6 - pkin(3) * t62 + t44 * t20 + t45 * t21 + t128 + (-t161 * g(2) - (t125 * t75 + t134 * t127) * t160 - pkin(3) * t57 + t134 * qJD(4) + (t152 - t155) * qJ(4)) * m(5) + (-t141 * g(2) - t100 * t42 + t175 * t14 + t176 * t15 - t56 * t173 + t2 * t45 + t3 * t44) * m(6); -t147 * t120 ^ 2 * mrSges(5,3) - t58 * t46 + t59 * t47 + t6 + t62 + (t14 * t59 - t15 * t58 + t138 + t42) * m(6) + (-t134 * t120 + t138 + t57) * m(5); Ifges(6,5) * t28 + Ifges(6,6) * t29 + Ifges(6,3) * qJDD(5) - t2 * mrSges(6,2) + t3 * mrSges(6,1) - t56 * (mrSges(6,1) * t59 + mrSges(6,2) * t58) - t59 * (Ifges(6,1) * t58 - t163) / 0.2e1 + t22 * t164 - qJD(5) * (Ifges(6,5) * t58 - Ifges(6,6) * t59) / 0.2e1 - t14 * t46 + t15 * t47 - g(3) * t172 + (t14 * t58 + t15 * t59) * mrSges(6,3) - (-Ifges(6,2) * t59 + t23 + t53) * t58 / 0.2e1 + (g(1) * t98 + g(2) * t97) * (mrSges(6,1) * t108 + mrSges(6,2) * t110);];
tau = t1;
