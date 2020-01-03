% Calculate vector of inverse dynamics joint torques for
% S4PRRR4
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
% Datum: 2019-12-31 16:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:32:30
% EndTime: 2019-12-31 16:32:34
% DurationCPUTime: 1.66s
% Computational Cost: add. (1034->220), mult. (2283->318), div. (0->0), fcn. (1374->8), ass. (0->104)
t79 = sin(qJ(3));
t81 = cos(qJ(3));
t57 = -mrSges(4,1) * t81 + mrSges(4,2) * t79;
t77 = qJ(3) + qJ(4);
t72 = sin(t77);
t73 = cos(t77);
t97 = t73 * mrSges(5,1) - t72 * mrSges(5,2);
t134 = t57 - t97;
t103 = qJD(2) * qJD(3);
t51 = qJDD(2) * t81 - t79 * t103;
t75 = pkin(7) + qJ(2);
t67 = sin(t75);
t68 = cos(t75);
t132 = g(1) * t68 + g(2) * t67;
t127 = t79 / 0.2e1;
t106 = qJDD(2) * pkin(2);
t104 = qJD(1) * qJD(3);
t32 = pkin(5) * t51 + t79 * qJDD(1) + t81 * t104;
t52 = qJDD(2) * t79 + t103 * t81;
t33 = -pkin(5) * t52 + qJDD(1) * t81 - t104 * t79;
t133 = t32 * t81 - t33 * t79;
t63 = pkin(3) * t81 + pkin(2);
t131 = m(4) * pkin(2) + m(5) * t63 + mrSges(3,1) - t134;
t82 = -pkin(6) - pkin(5);
t130 = -m(4) * pkin(5) + m(5) * t82 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t78 = sin(qJ(4));
t80 = cos(qJ(4));
t48 = t78 * t81 + t79 * t80;
t46 = t48 * qJD(2);
t128 = t46 / 0.2e1;
t125 = m(2) + m(3);
t122 = Ifges(4,4) * t79;
t121 = Ifges(4,4) * t81;
t120 = Ifges(4,2) * t81;
t111 = qJD(1) * t79;
t60 = t82 * t81;
t42 = -qJD(2) * t60 + t111;
t116 = t42 * t78;
t59 = t82 * t79;
t71 = t81 * qJD(1);
t41 = qJD(2) * t59 + t71;
t36 = qJD(3) * pkin(3) + t41;
t17 = t36 * t80 - t116;
t119 = t17 * mrSges(5,3);
t115 = t42 * t80;
t114 = t46 * mrSges(5,3);
t113 = t46 * Ifges(5,4);
t110 = qJD(2) * t79;
t109 = qJD(2) * t81;
t108 = qJD(3) * t79;
t107 = qJD(3) * t81;
t102 = pkin(3) * t108;
t99 = qJD(3) * t82;
t95 = mrSges(4,1) * t79 + mrSges(4,2) * t81;
t94 = mrSges(5,1) * t72 + mrSges(5,2) * t73;
t93 = t120 + t122;
t92 = Ifges(4,5) * t81 - Ifges(4,6) * t79;
t18 = t36 * t78 + t115;
t53 = -pkin(5) * t110 + t71;
t54 = pkin(5) * t109 + t111;
t91 = t53 * t81 + t54 * t79;
t30 = t59 * t80 + t60 * t78;
t31 = t59 * t78 - t60 * t80;
t47 = -t78 * t79 + t80 * t81;
t90 = pkin(2) * t95;
t89 = t79 * (Ifges(4,1) * t81 - t122);
t88 = t48 * qJD(4);
t87 = t47 * qJD(4);
t15 = qJD(2) * t87 + t51 * t78 + t52 * t80;
t16 = -qJD(2) * t88 + t51 * t80 - t52 * t78;
t19 = qJDD(3) * pkin(3) - pkin(6) * t52 + t33;
t25 = pkin(6) * t51 + t32;
t2 = qJD(4) * t17 + t19 * t78 + t25 * t80;
t45 = t47 * qJD(2);
t76 = qJD(3) + qJD(4);
t23 = t45 * Ifges(5,2) + t76 * Ifges(5,6) + t113;
t40 = Ifges(5,4) * t45;
t24 = t46 * Ifges(5,1) + t76 * Ifges(5,5) + t40;
t3 = -qJD(4) * t18 + t19 * t80 - t25 * t78;
t58 = t63 * qJD(2);
t74 = qJDD(3) + qJDD(4);
t83 = t3 * mrSges(5,1) - t2 * mrSges(5,2) + t45 * t119 + t23 * t128 + t58 * (mrSges(5,1) * t46 + mrSges(5,2) * t45) - t46 * (Ifges(5,1) * t45 - t113) / 0.2e1 + Ifges(5,6) * t16 + Ifges(5,5) * t15 - t76 * (Ifges(5,5) * t45 - Ifges(5,6) * t46) / 0.2e1 + Ifges(5,3) * t74 - (-Ifges(5,2) * t46 + t24 + t40) * t45 / 0.2e1;
t65 = Ifges(4,4) * t109;
t56 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t109;
t55 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t110;
t50 = t81 * t99;
t49 = t79 * t99;
t44 = Ifges(4,1) * t110 + Ifges(4,5) * qJD(3) + t65;
t43 = Ifges(4,6) * qJD(3) + qJD(2) * t93;
t39 = -pkin(3) * t51 - t106;
t38 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t52;
t37 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t51;
t35 = mrSges(5,1) * t76 - t114;
t34 = -mrSges(5,2) * t76 + mrSges(5,3) * t45;
t29 = -qJD(3) * t48 - t88;
t28 = qJD(3) * t47 + t87;
t27 = -mrSges(5,1) * t45 + mrSges(5,2) * t46;
t21 = t41 * t80 - t116;
t20 = -t41 * t78 - t115;
t10 = -mrSges(5,2) * t74 + mrSges(5,3) * t16;
t9 = mrSges(5,1) * t74 - mrSges(5,3) * t15;
t6 = -qJD(4) * t31 - t49 * t78 + t50 * t80;
t5 = qJD(4) * t30 + t49 * t80 + t50 * t78;
t1 = [t48 * t10 + t28 * t34 + t29 * t35 + t79 * t37 + t81 * t38 + t47 * t9 + t125 * qJDD(1) + (-t55 * t79 + t56 * t81) * qJD(3) + m(4) * (t32 * t79 + t33 * t81 + (-t53 * t79 + t54 * t81) * qJD(3)) + m(5) * (t17 * t29 + t18 * t28 + t2 * t48 + t3 * t47) + (-m(4) - m(5) - t125) * g(3); t81 * (Ifges(4,4) * t52 + Ifges(4,2) * t51) / 0.2e1 + t27 * t102 + t52 * (Ifges(4,1) * t79 + t121) / 0.2e1 - t28 * t119 - t57 * t106 + t44 * t107 / 0.2e1 - t43 * t108 / 0.2e1 - t90 * t103 + (Ifges(5,1) * t28 + Ifges(5,4) * t29) * t128 + (Ifges(4,1) * t52 + Ifges(4,4) * t51) * t127 + qJD(3) ^ 2 * t92 / 0.2e1 + t51 * t93 / 0.2e1 + (t81 * (-Ifges(4,2) * t79 + t121) + t89) * t103 / 0.2e1 + m(5) * (-t102 * t58 + t17 * t6 + t18 * t5 + t2 * t31 + t3 * t30 - t39 * t63) + (0.2e1 * Ifges(4,5) * t127 + Ifges(4,6) * t81) * qJDD(3) + (-t107 * t53 - t108 * t54 + t133) * mrSges(4,3) + (-t55 * t107 - t56 * t108 + m(4) * (-qJD(3) * t91 + t133) - t79 * t38 + t81 * t37) * pkin(5) + (m(4) * t106 + t51 * mrSges(4,1) - t52 * mrSges(4,2)) * pkin(2) + (t130 * t68 + t131 * t67) * g(1) + (t130 * t67 - t131 * t68) * g(2) + (mrSges(5,2) * t39 - mrSges(5,3) * t3 + Ifges(5,1) * t15 + Ifges(5,4) * t16 + Ifges(5,5) * t74) * t48 + (-mrSges(5,1) * t39 + mrSges(5,3) * t2 + Ifges(5,4) * t15 + Ifges(5,2) * t16 + Ifges(5,6) * t74) * t47 + t76 * (Ifges(5,5) * t28 + Ifges(5,6) * t29) / 0.2e1 - t63 * (-mrSges(5,1) * t16 + mrSges(5,2) * t15) - t58 * (-mrSges(5,1) * t29 + mrSges(5,2) * t28) + t45 * (Ifges(5,4) * t28 + Ifges(5,2) * t29) / 0.2e1 + t30 * t9 + t31 * t10 + t5 * t34 + t6 * t35 + t28 * t24 / 0.2e1 + t29 * t23 / 0.2e1 + t18 * t29 * mrSges(5,3) + Ifges(3,3) * qJDD(2); t134 * g(3) - m(5) * (t17 * t20 + t18 * t21) + (-qJD(3) * t92 / 0.2e1 + t43 * t127 + (t120 * t127 - t89 / 0.2e1 + t90) * qJD(2) + (m(5) * t58 - t27) * t79 * pkin(3) + t91 * mrSges(4,3) - (t44 + t65) * t81 / 0.2e1) * qJD(2) + (t78 * t10 + t80 * t9 + (-g(3) * t81 + t132 * t79 + t2 * t78 + t3 * t80) * m(5) + (t80 * t34 - t78 * t35 + (-t17 * t78 + t18 * t80) * m(5)) * qJD(4)) * pkin(3) + t83 + Ifges(4,6) * t51 + Ifges(4,5) * t52 + t54 * t55 - t53 * t56 - t32 * mrSges(4,2) + t33 * mrSges(4,1) - t21 * t34 - t20 * t35 + t18 * t114 + Ifges(4,3) * qJDD(3) + t132 * (t94 + t95); (t35 + t114) * t18 + t83 - g(3) * t97 - t17 * t34 + t132 * t94;];
tau = t1;
