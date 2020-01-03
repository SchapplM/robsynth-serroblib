% Calculate vector of inverse dynamics joint torques for
% S4RPRR4
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
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:14
% EndTime: 2019-12-31 16:50:21
% DurationCPUTime: 3.38s
% Computational Cost: add. (1297->300), mult. (2798->440), div. (0->0), fcn. (1564->10), ass. (0->142)
t81 = cos(qJ(4));
t133 = qJD(3) * t81;
t79 = sin(qJ(3));
t137 = qJD(1) * t79;
t78 = sin(qJ(4));
t50 = -t137 * t78 + t133;
t184 = -t50 / 0.2e1;
t135 = qJD(3) * t78;
t51 = t137 * t81 + t135;
t183 = -t51 / 0.2e1;
t82 = cos(qJ(3));
t128 = t82 * qJD(1);
t63 = qJD(4) - t128;
t182 = -t63 / 0.2e1;
t75 = qJ(1) + pkin(7);
t71 = sin(t75);
t72 = cos(t75);
t181 = g(1) * t72 + g(2) * t71;
t179 = -m(5) - m(4);
t127 = qJD(1) * qJD(3);
t55 = qJDD(1) * t79 + t127 * t82;
t24 = qJD(4) * t50 + qJDD(3) * t78 + t55 * t81;
t25 = -qJD(4) * t51 + qJDD(3) * t81 - t55 * t78;
t5 = -mrSges(5,1) * t25 + mrSges(5,2) * t24;
t178 = -qJDD(3) * mrSges(4,1) + mrSges(4,3) * t55 + t5;
t155 = Ifges(4,4) * t79;
t177 = -qJD(1) / 0.2e1;
t119 = mrSges(4,3) * t137;
t176 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t50 + mrSges(5,2) * t51 + t119;
t47 = Ifges(5,4) * t50;
t18 = Ifges(5,1) * t51 + Ifges(5,5) * t63 + t47;
t70 = Ifges(4,4) * t128;
t175 = Ifges(4,1) * t137 + Ifges(4,5) * qJD(3) + t81 * t18 + t70;
t76 = sin(pkin(7));
t67 = pkin(1) * t76 + pkin(5);
t58 = t67 * qJD(1);
t39 = qJD(2) * t79 + t58 * t82;
t136 = qJD(3) * t39;
t56 = t67 * qJDD(1);
t15 = qJDD(2) * t82 - t56 * t79 - t136;
t32 = qJD(3) * pkin(6) + t39;
t77 = cos(pkin(7));
t162 = pkin(1) * t77;
t111 = pkin(3) * t82 + pkin(6) * t79;
t95 = -pkin(2) - t111;
t46 = t95 - t162;
t33 = t46 * qJD(1);
t10 = -t32 * t78 + t33 * t81;
t11 = t32 * t81 + t33 * t78;
t129 = qJD(4) * t81;
t131 = qJD(4) * t78;
t174 = -t10 * t129 - t11 * t131;
t132 = qJD(3) * t82;
t134 = qJD(3) * t79;
t14 = qJD(2) * t132 + t79 * qJDD(2) - t134 * t58 + t82 * t56;
t173 = t14 * t82 - t15 * t79;
t54 = qJDD(1) * t82 - t127 * t79;
t48 = qJDD(4) - t54;
t8 = mrSges(5,1) * t48 - mrSges(5,3) * t24;
t9 = -mrSges(5,2) * t48 + mrSges(5,3) * t25;
t172 = -t78 * t8 + t81 * t9;
t102 = t82 * Ifges(4,2) + t155;
t170 = Ifges(4,6) * qJD(3) / 0.2e1 + qJD(1) * t102 / 0.2e1 + Ifges(5,5) * t183 + Ifges(5,6) * t184 + Ifges(5,3) * t182;
t108 = t82 * mrSges(4,1) - mrSges(4,2) * t79;
t146 = t79 * mrSges(5,3);
t169 = -mrSges(3,1) - t108 - t146;
t168 = t24 / 0.2e1;
t167 = t25 / 0.2e1;
t166 = t48 / 0.2e1;
t164 = t51 / 0.2e1;
t163 = t81 / 0.2e1;
t80 = sin(qJ(1));
t161 = pkin(1) * t80;
t12 = qJDD(3) * pkin(6) + t14;
t68 = -pkin(2) - t162;
t57 = t68 * qJDD(1);
t23 = -pkin(3) * t54 - pkin(6) * t55 + t57;
t2 = -qJD(4) * t11 - t12 * t78 + t23 * t81;
t158 = t2 * t78;
t83 = cos(qJ(1));
t74 = t83 * pkin(1);
t154 = Ifges(5,4) * t51;
t153 = Ifges(5,4) * t78;
t152 = Ifges(5,4) * t81;
t17 = Ifges(5,2) * t50 + Ifges(5,6) * t63 + t154;
t148 = t78 * t17;
t147 = t78 * t82;
t145 = t81 * mrSges(5,3);
t143 = t81 * t82;
t142 = t82 * mrSges(5,3);
t141 = t82 * Ifges(4,4);
t130 = qJD(4) * t79;
t126 = Ifges(5,5) * t24 + Ifges(5,6) * t25 + Ifges(5,3) * t48;
t124 = mrSges(4,3) * t128;
t123 = t67 * t134;
t118 = -t148 / 0.2e1;
t117 = m(5) * pkin(6) + mrSges(5,3);
t110 = pkin(3) * t79 - pkin(6) * t82;
t1 = qJD(4) * t10 + t12 * t81 + t23 * t78;
t109 = t1 * t81 - t158;
t107 = mrSges(4,1) * t79 + mrSges(4,2) * t82;
t106 = -t81 * mrSges(5,1) + t78 * mrSges(5,2);
t105 = mrSges(5,1) * t78 + mrSges(5,2) * t81;
t104 = Ifges(5,1) * t81 - t153;
t103 = Ifges(5,1) * t78 + t152;
t101 = -Ifges(5,2) * t78 + t152;
t100 = Ifges(5,2) * t81 + t153;
t99 = Ifges(4,5) * t82 - Ifges(4,6) * t79;
t98 = Ifges(5,5) * t81 - Ifges(5,6) * t78;
t97 = Ifges(5,5) * t78 + Ifges(5,6) * t81;
t38 = qJD(2) * t82 - t58 * t79;
t27 = t143 * t67 + t46 * t78;
t26 = -t147 * t67 + t46 * t81;
t31 = -qJD(3) * pkin(3) - t38;
t94 = t31 * t105;
t93 = t68 * qJD(1) * t107;
t92 = t79 * (Ifges(4,1) * t82 - t155);
t91 = m(5) * pkin(3) - t106;
t90 = -t130 * t78 + t132 * t81;
t89 = t129 * t79 + t132 * t78;
t88 = Ifges(5,5) * t79 + t104 * t82;
t87 = Ifges(5,6) * t79 + t101 * t82;
t86 = Ifges(5,3) * t79 + t82 * t98;
t61 = -qJD(3) * mrSges(4,2) + t124;
t53 = t110 * qJD(3);
t52 = t110 * qJD(1);
t45 = t105 * t79;
t40 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t54;
t37 = t143 * t72 + t71 * t78;
t36 = -t147 * t72 + t71 * t81;
t35 = -t143 * t71 + t72 * t78;
t34 = t147 * t71 + t72 * t81;
t30 = mrSges(5,1) * t63 - mrSges(5,3) * t51;
t29 = -mrSges(5,2) * t63 + mrSges(5,3) * t50;
t20 = t38 * t81 + t52 * t78;
t19 = -t38 * t78 + t52 * t81;
t13 = -qJDD(3) * pkin(3) - t15;
t7 = -qJD(4) * t27 + t123 * t78 + t53 * t81;
t6 = qJD(4) * t26 - t123 * t81 + t53 * t78;
t4 = Ifges(5,1) * t24 + Ifges(5,4) * t25 + Ifges(5,5) * t48;
t3 = Ifges(5,4) * t24 + Ifges(5,2) * t25 + Ifges(5,6) * t48;
t16 = [t82 * (Ifges(4,4) * t55 + Ifges(4,2) * t54) / 0.2e1 + t55 * t79 * Ifges(4,1) + (-t134 * t39 + t173 - t181) * mrSges(4,3) + (-m(3) * t74 - mrSges(2,1) * t83 - t37 * mrSges(5,1) + mrSges(2,2) * t80 + mrSges(3,2) * t71 - t36 * mrSges(5,2) + t179 * (t72 * pkin(2) + t71 * pkin(5) + t74) + (-m(5) * t111 + t169) * t72) * g(2) + (m(3) * t161 + mrSges(2,1) * t80 - t35 * mrSges(5,1) + mrSges(2,2) * t83 + mrSges(3,2) * t72 - t34 * mrSges(5,2) + t179 * (t72 * pkin(5) - t161) + (m(4) * pkin(2) - m(5) * t95 - t169) * t71) * g(1) + (t10 * mrSges(5,1) - t11 * mrSges(5,2) - t170) * t134 + m(5) * (t1 * t27 + t10 * t7 + t11 * t6 + t2 * t26) + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t77 * mrSges(3,1) - 0.2e1 * t76 * mrSges(3,2) + m(3) * (t76 ^ 2 + t77 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t155 + t102) * t54 / 0.2e1 + (m(4) * t68 - t108) * t57 + qJDD(3) * (Ifges(4,5) * t79 + t82 * Ifges(4,6)) + t68 * (-mrSges(4,1) * t54 + mrSges(4,2) * t55) - t78 * t79 * t3 / 0.2e1 + (t82 * (-Ifges(4,2) * t79 + t141) + t92) * t127 / 0.2e1 - (t81 * t17 + t78 * t18) * t130 / 0.2e1 + (-t38 * mrSges(4,3) + t175 / 0.2e1 + t118) * t132 + t55 * t141 / 0.2e1 + t13 * t45 + t26 * t8 + t27 * t9 + t6 * t29 + t7 * t30 + (m(4) * ((-t38 * t82 - t39 * t79) * qJD(3) + t173) + t178 * t79 + t176 * t132 + t82 * t40 + m(5) * (t13 * t79 + t132 * t31)) * t67 + (-t10 * t90 - t11 * t89) * mrSges(5,3) + (qJD(3) * t88 - t103 * t130) * t164 + (-Ifges(5,3) * t82 + t79 * t98) * t166 + (-Ifges(5,6) * t82 + t101 * t79) * t167 + (-Ifges(5,5) * t82 + t104 * t79) * t168 + t31 * (mrSges(5,1) * t89 + mrSges(5,2) * t90) + qJD(3) ^ 2 * t99 / 0.2e1 + t79 * t4 * t163 + qJD(3) * t93 - t61 * t123 - t82 * t126 / 0.2e1 + t50 * (qJD(3) * t87 - t100 * t130) / 0.2e1 + t63 * (qJD(3) * t86 - t130 * t97) / 0.2e1 + t2 * (-mrSges(5,1) * t82 - t145 * t79) + t1 * (t82 * mrSges(5,2) - t146 * t78); m(3) * qJDD(2) + (-m(3) + t179) * g(3) + ((t81 * t29 - t78 * t30 + t61) * qJD(3) + m(4) * (t15 + t136) + m(5) * (-t10 * t135 + t11 * t133 - t13) - t178) * t82 + (t40 + (-t78 * t29 - t81 * t30) * qJD(4) + t176 * qJD(3) + m(4) * (-qJD(3) * t38 + t14) + m(5) * (qJD(3) * t31 + t109 + t174) + t172) * t79; (t148 / 0.2e1 - t94) * t128 + (t118 + t94) * qJD(4) + (-t61 + t124) * t38 + (-t117 * t82 + t79 * t91 + t107) * t181 + (-m(5) * t31 + t119 - t176) * t39 + (-t10 * (mrSges(5,1) * t79 - t142 * t81) - t11 * (-mrSges(5,2) * t79 - t142 * t78) - t93 + t92 * t177) * qJD(1) + (-t158 + t174) * mrSges(5,3) - (-Ifges(4,2) * t137 + t175 + t70) * t128 / 0.2e1 + t170 * t137 + (-t30 * t129 - t29 * t131 + m(5) * ((-t10 * t81 - t11 * t78) * qJD(4) + t109) + t172) * pkin(6) + (-pkin(3) * t13 - t10 * t19 - t11 * t20) * m(5) + (-t117 * t79 - t82 * t91 - t108) * g(3) + t78 * t4 / 0.2e1 + Ifges(4,6) * t54 + Ifges(4,5) * t55 + (t101 * t50 + t104 * t51 + t63 * t98) * qJD(4) / 0.2e1 + t1 * t145 - t20 * t29 - t19 * t30 - t14 * mrSges(4,2) + t15 * mrSges(4,1) - pkin(3) * t5 + t3 * t163 + t97 * t166 + t100 * t167 + t103 * t168 + t13 * t106 + (t50 * t87 + t51 * t88 + t63 * t86) * t177 + Ifges(4,3) * qJDD(3) - t99 * t127 / 0.2e1 + t18 * t129 / 0.2e1; -t1 * mrSges(5,2) + t2 * mrSges(5,1) - t31 * (mrSges(5,1) * t51 + mrSges(5,2) * t50) + (Ifges(5,1) * t50 - t154) * t183 + t17 * t164 + (Ifges(5,5) * t50 - Ifges(5,6) * t51) * t182 - t10 * t29 + t11 * t30 - g(1) * (mrSges(5,1) * t36 - mrSges(5,2) * t37) - g(2) * (-mrSges(5,1) * t34 + mrSges(5,2) * t35) + g(3) * t45 + (t10 * t50 + t11 * t51) * mrSges(5,3) + t126 + (-Ifges(5,2) * t51 + t18 + t47) * t184;];
tau = t16;
