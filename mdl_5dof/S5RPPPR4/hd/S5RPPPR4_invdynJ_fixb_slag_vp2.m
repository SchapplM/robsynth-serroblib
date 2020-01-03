% Calculate vector of inverse dynamics joint torques for
% S5RPPPR4
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
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
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
% Datum: 2019-12-31 17:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:45:06
% EndTime: 2019-12-31 17:45:09
% DurationCPUTime: 1.92s
% Computational Cost: add. (1076->204), mult. (1949->268), div. (0->0), fcn. (1134->12), ass. (0->93)
t79 = sin(pkin(8));
t81 = cos(pkin(8));
t84 = sin(qJ(5));
t86 = cos(qJ(5));
t135 = -t79 * t84 + t81 * t86;
t39 = t135 * qJD(5);
t80 = sin(pkin(7));
t56 = pkin(1) * t80 + qJ(3);
t136 = t79 ^ 2 + t81 ^ 2;
t104 = qJD(1) * qJD(4);
t82 = cos(pkin(7));
t62 = -pkin(1) * t82 - pkin(2);
t53 = -qJ(4) + t62;
t133 = qJDD(1) * t53;
t30 = qJDD(3) - t104 + t133;
t21 = -qJDD(2) * t79 + t81 * t30;
t22 = t81 * qJDD(2) + t79 * t30;
t93 = t21 * t81 + t22 * t79;
t75 = qJ(1) + pkin(7);
t68 = sin(t75);
t70 = cos(t75);
t134 = -g(1) * t68 + g(2) * t70;
t105 = m(4) + m(5) + m(6);
t100 = t136 * mrSges(5,3);
t132 = mrSges(3,1) + mrSges(5,3) + mrSges(6,3) - mrSges(4,2);
t124 = pkin(4) * t79;
t74 = pkin(8) + qJ(5);
t67 = sin(t74);
t69 = cos(t74);
t94 = -t67 * mrSges(6,1) - t69 * mrSges(6,2);
t96 = mrSges(5,1) * t79 + mrSges(5,2) * t81;
t131 = -m(6) * t124 + mrSges(3,2) - mrSges(4,3) + t94 - t96;
t106 = qJDD(1) * t81;
t17 = -pkin(6) * t106 + t21;
t107 = qJDD(1) * t79;
t18 = -pkin(6) * t107 + t22;
t108 = qJD(1) * t81;
t42 = qJD(1) * t53 + qJD(3);
t23 = -qJD(2) * t79 + t81 * t42;
t19 = -pkin(6) * t108 + t23;
t109 = qJD(1) * t79;
t24 = t81 * qJD(2) + t79 * t42;
t20 = -pkin(6) * t109 + t24;
t3 = t19 * t86 - t20 * t84;
t1 = qJD(5) * t3 + t17 * t84 + t18 * t86;
t4 = t19 * t84 + t20 * t86;
t2 = -qJD(5) * t4 + t17 * t86 - t18 * t84;
t44 = -t86 * t79 - t84 * t81;
t40 = t44 * qJD(5);
t130 = t1 * t44 - t135 * t2 - t3 * t40 - t4 * t39;
t50 = t56 * qJD(1);
t48 = qJD(4) + t50;
t34 = pkin(4) * t109 + t48;
t37 = t44 * qJD(1);
t38 = t108 * t86 - t109 * t84;
t129 = -m(4) * t50 - m(5) * t48 - m(6) * t34 + mrSges(6,1) * t37 - mrSges(6,2) * t38 - t96 * qJD(1);
t127 = t38 / 0.2e1;
t85 = sin(qJ(1));
t125 = pkin(1) * t85;
t87 = cos(qJ(1));
t71 = t87 * pkin(1);
t117 = -pkin(6) + t53;
t116 = Ifges(6,4) * t38;
t111 = mrSges(5,1) * t107 + mrSges(5,2) * t106;
t110 = pkin(1) * qJDD(1);
t77 = qJD(1) * qJD(3);
t63 = t80 * t110;
t43 = -qJDD(1) * qJ(3) - t63 - t77;
t15 = qJD(1) * t40 + qJDD(1) * t135;
t16 = -qJD(1) * t39 + qJDD(1) * t44;
t98 = -t16 * mrSges(6,1) + t15 * mrSges(6,2);
t41 = qJDD(4) - t43;
t97 = -g(1) * t70 - g(2) * t68;
t92 = -t23 * t81 - t24 * t79;
t35 = t117 * t79;
t36 = t117 * t81;
t10 = t35 * t86 + t36 * t84;
t9 = -t35 * t84 + t36 * t86;
t88 = qJD(1) ^ 2;
t83 = -pkin(6) - qJ(4);
t49 = t56 + t124;
t47 = qJDD(1) * t62 + qJDD(3);
t33 = Ifges(6,4) * t37;
t29 = pkin(4) * t107 + t41;
t28 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t38;
t27 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t37;
t12 = Ifges(6,1) * t38 + Ifges(6,5) * qJD(5) + t33;
t11 = Ifges(6,2) * t37 + Ifges(6,6) * qJD(5) + t116;
t8 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t16;
t7 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t15;
t6 = -qJD(4) * t135 - qJD(5) * t10;
t5 = qJD(4) * t44 + qJD(5) * t9;
t13 = [t56 * t111 - (Ifges(5,4) * t81 - Ifges(5,2) * t79) * t107 - 0.2e1 * mrSges(3,2) * t63 + t49 * t98 + t41 * t96 + 0.2e1 * t82 * mrSges(3,1) * t110 + (m(3) * (t80 ^ 2 + t82 ^ 2) * pkin(1) ^ 2 + t56 * mrSges(4,3) + t62 * mrSges(4,2) + Ifges(4,1) + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + (-t43 + t77) * mrSges(4,3) + (-mrSges(6,1) * t29 + Ifges(6,4) * t15 + Ifges(6,2) * t16 + Ifges(6,6) * qJDD(5)) * t44 - t100 * t133 + m(6) * (t1 * t10 + t2 * t9 + t29 * t49 + t3 * t6 + t4 * t5) + m(4) * (-t43 * t56 + t47 * t62) + (t104 * t136 - t93) * mrSges(5,3) + (m(3) * t125 + t85 * mrSges(2,1) + t87 * mrSges(2,2) - t105 * (t70 * qJ(3) - t125) + t131 * t70 + (-m(5) * (-pkin(2) - qJ(4)) - m(6) * (-pkin(2) + t83) + m(4) * pkin(2) + t132) * t68) * g(1) + (-m(3) * t71 - t87 * mrSges(2,1) + t85 * mrSges(2,2) - t105 * (t70 * pkin(2) + t68 * qJ(3) + t71) + (-m(5) * qJ(4) + m(6) * t83 - t132) * t70 + t131 * t68) * g(2) + t130 * mrSges(6,3) - t129 * qJD(3) + (mrSges(6,2) * t29 + Ifges(6,1) * t15 + Ifges(6,4) * t16 + Ifges(6,5) * qJDD(5)) * t135 + t47 * mrSges(4,2) - t39 * t11 / 0.2e1 + t37 * (Ifges(6,4) * t40 - Ifges(6,2) * t39) / 0.2e1 + t34 * (t39 * mrSges(6,1) + t40 * mrSges(6,2)) + qJD(5) * (Ifges(6,5) * t40 - Ifges(6,6) * t39) / 0.2e1 + t40 * t12 / 0.2e1 + t5 * t27 + t6 * t28 + t9 * t7 + t10 * t8 + m(5) * (t92 * qJD(4) + t41 * t56 + t93 * t53) + (Ifges(5,1) * t81 - Ifges(5,4) * t79) * t106 + (Ifges(6,1) * t40 - Ifges(6,4) * t39) * t127; t40 * t27 - t39 * t28 + t44 * t7 + t135 * t8 + (m(3) + m(4)) * qJDD(2) + m(5) * (-t21 * t79 + t22 * t81) + m(6) * (t1 * t135 + t2 * t44 - t3 * t39 + t4 * t40) + (-m(3) - t105) * g(3); -t88 * mrSges(4,3) + t39 * t27 + t40 * t28 - t44 * t8 + t135 * t7 + (mrSges(4,2) - t100) * qJDD(1) + m(4) * t47 + m(5) * t93 - m(6) * t130 + t129 * qJD(1) + t134 * t105; -t88 * t100 - t37 * t27 + t38 * t28 + t111 + t98 + (t3 * t38 - t4 * t37 + t29 + t97) * m(6) + (-qJD(1) * t92 + t41 + t97) * m(5); Ifges(6,5) * t15 + Ifges(6,6) * t16 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t34 * (mrSges(6,1) * t38 + mrSges(6,2) * t37) - t38 * (Ifges(6,1) * t37 - t116) / 0.2e1 + t11 * t127 - qJD(5) * (Ifges(6,5) * t37 - Ifges(6,6) * t38) / 0.2e1 - t3 * t27 + t4 * t28 - g(3) * t94 + (t3 * t37 + t38 * t4) * mrSges(6,3) + t134 * (mrSges(6,1) * t69 - mrSges(6,2) * t67) - (-Ifges(6,2) * t38 + t12 + t33) * t37 / 0.2e1;];
tau = t13;
