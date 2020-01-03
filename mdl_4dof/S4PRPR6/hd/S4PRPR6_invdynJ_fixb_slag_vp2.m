% Calculate vector of inverse dynamics joint torques for
% S4PRPR6
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
%   pkin=[a2,a3,a4,d2,d4,theta1,theta3]';
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
% Datum: 2019-12-31 16:24
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:24:20
% EndTime: 2019-12-31 16:24:24
% DurationCPUTime: 1.94s
% Computational Cost: add. (705->180), mult. (1721->260), div. (0->0), fcn. (1107->10), ass. (0->87)
t59 = sin(pkin(7));
t61 = cos(pkin(7));
t98 = t59 ^ 2 + t61 ^ 2;
t120 = t98 * mrSges(4,3);
t119 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t67 = cos(qJ(2));
t81 = -qJD(1) * t67 + qJD(3);
t50 = pkin(3) * t61 + pkin(2);
t58 = pkin(7) + qJ(4);
t53 = sin(t58);
t54 = cos(t58);
t77 = -mrSges(4,1) * t61 + mrSges(4,2) * t59;
t118 = m(4) * pkin(2) + m(5) * t50 + mrSges(5,1) * t54 - mrSges(5,2) * t53 + mrSges(3,1) - t77;
t99 = pkin(5) + qJ(3);
t117 = -m(4) * qJ(3) - m(5) * t99 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t60 = sin(pkin(6));
t62 = cos(pkin(6));
t116 = g(1) * t62 + g(2) * t60;
t45 = t99 * t59;
t46 = t99 * t61;
t64 = sin(qJ(4));
t66 = cos(qJ(4));
t21 = -t45 * t64 + t46 * t66;
t41 = t59 * t66 + t61 * t64;
t74 = t41 * t67;
t115 = qJD(1) * t74 - qJD(3) * t41 - qJD(4) * t21;
t20 = -t45 * t66 - t46 * t64;
t76 = t59 * t64 - t61 * t66;
t114 = qJD(4) * t20 - t81 * t76;
t36 = t76 * qJD(4);
t16 = -qJD(2) * t36 + qJDD(2) * t41;
t37 = t41 * qJD(4);
t17 = -qJD(2) * t37 - qJDD(2) * t76;
t3 = -t17 * mrSges(5,1) + t16 * mrSges(5,2);
t94 = qJDD(2) * t61;
t95 = qJDD(2) * t59;
t39 = -mrSges(4,1) * t94 + mrSges(4,2) * t95;
t113 = t3 + t39;
t65 = sin(qJ(2));
t112 = t65 * t116;
t34 = t76 * qJD(2);
t35 = t41 * qJD(2);
t111 = mrSges(5,1) * t34 + mrSges(5,2) * t35 + t77 * qJD(2);
t107 = t35 / 0.2e1;
t106 = m(4) + m(5);
t103 = Ifges(5,4) * t35;
t101 = t67 * t53;
t100 = t67 * t54;
t93 = qJD(1) * qJD(2);
t52 = t67 * t93;
t44 = t65 * qJDD(1) + t52;
t97 = qJD(1) * t65;
t48 = qJD(2) * qJ(3) + t97;
t89 = t98 * t48;
t88 = t98 * t65;
t87 = t98 * t67;
t51 = t65 * t93;
t83 = pkin(5) * qJD(2) + t48;
t29 = t44 + t119;
t82 = pkin(5) * qJDD(2) + t29;
t43 = qJDD(1) * t67 - t51;
t27 = t83 * t59;
t28 = t83 * t61;
t6 = -t27 * t66 - t28 * t64;
t7 = -t27 * t64 + t28 * t66;
t75 = qJDD(3) - t43;
t71 = (-qJD(2) * pkin(2) + t81) * t65 + t48 * t87;
t68 = qJD(2) ^ 2;
t38 = -qJD(2) * t50 + t81;
t33 = -qJDD(2) * pkin(2) + t75;
t32 = Ifges(5,4) * t34;
t31 = t76 * t65;
t30 = t41 * t65;
t24 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t35;
t23 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t34;
t22 = -qJDD(2) * t50 + t75;
t19 = t82 * t61;
t18 = t82 * t59;
t13 = Ifges(5,1) * t35 + Ifges(5,5) * qJD(4) - t32;
t12 = -Ifges(5,2) * t34 + Ifges(5,6) * qJD(4) + t103;
t11 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t17;
t10 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t16;
t9 = -qJD(2) * t74 + t65 * t36;
t8 = -t34 * t67 - t37 * t65;
t2 = -qJD(4) * t7 - t18 * t66 - t19 * t64;
t1 = qJD(4) * t6 - t18 * t64 + t19 * t66;
t4 = [m(2) * qJDD(1) - t30 * t10 - t31 * t11 + t8 * t23 + t9 * t24 + (qJDD(2) * mrSges(3,1) - t68 * mrSges(3,2) - t113) * t67 + (-t68 * mrSges(3,1) - qJDD(2) * mrSges(3,2) + t111 * qJD(2)) * t65 + m(3) * (t43 * t67 + t44 * t65) + m(4) * (qJD(2) * t71 + t29 * t88 - t33 * t67) + m(5) * (qJD(2) * t38 * t65 - t1 * t31 - t2 * t30 - t22 * t67 + t6 * t9 + t7 * t8) + (-m(2) - m(3) - t106) * g(3) + (qJDD(2) * t88 + t68 * t87) * mrSges(4,3); (t6 * t36 - t7 * t37) * mrSges(5,3) + (-Ifges(5,1) * t36 - Ifges(5,4) * t37) * t107 + t38 * (mrSges(5,1) * t37 - mrSges(5,2) * t36) - t34 * (-Ifges(5,4) * t36 - Ifges(5,2) * t37) / 0.2e1 + qJD(4) * (-Ifges(5,5) * t36 - Ifges(5,6) * t37) / 0.2e1 - (-mrSges(5,1) * t22 + mrSges(5,3) * t1 + Ifges(5,4) * t16 + Ifges(5,2) * t17 + Ifges(5,6) * qJDD(4)) * t76 + (t29 - t52 + t119) * t120 + (qJ(3) * t29 * t98 - pkin(2) * t33 - t71 * qJD(1) + qJD(3) * t89) * m(4) + (Ifges(4,4) * t59 + Ifges(4,2) * t61) * t94 + (Ifges(4,1) * t59 + Ifges(4,4) * t61) * t95 + (t117 * g(3) + t116 * t118) * t65 + (-t118 * g(3) + t116 * t117) * t67 - t50 * t3 - pkin(2) * t39 - t36 * t13 / 0.2e1 - t37 * t12 / 0.2e1 + t20 * t10 + t21 * t11 + (t52 - t44) * mrSges(3,2) + (t51 + t43) * mrSges(3,1) + t33 * t77 + Ifges(3,3) * qJDD(2) - t111 * t97 + t114 * t23 + t115 * t24 + (t1 * t21 + t114 * t7 + t115 * t6 + t2 * t20 - t22 * t50 - t38 * t97) * m(5) + (mrSges(5,2) * t22 - mrSges(5,3) * t2 + Ifges(5,1) * t16 + Ifges(5,4) * t17 + Ifges(5,5) * qJDD(4)) * t41; -t68 * t120 + t106 * t67 * g(3) + t34 * t23 + t35 * t24 + (t34 * t7 + t35 * t6 - t112 + t22) * m(5) + (-qJD(2) * t89 - t112 + t33) * m(4) + t113; Ifges(5,5) * t16 + Ifges(5,6) * t17 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t38 * (mrSges(5,1) * t35 - mrSges(5,2) * t34) - t35 * (-Ifges(5,1) * t34 - t103) / 0.2e1 + t12 * t107 - qJD(4) * (-Ifges(5,5) * t34 - Ifges(5,6) * t35) / 0.2e1 - t6 * t23 + t7 * t24 - g(1) * ((-t101 * t62 + t54 * t60) * mrSges(5,1) + (-t100 * t62 - t53 * t60) * mrSges(5,2)) - g(2) * ((-t101 * t60 - t54 * t62) * mrSges(5,1) + (-t100 * t60 + t53 * t62) * mrSges(5,2)) - g(3) * (-mrSges(5,1) * t53 - mrSges(5,2) * t54) * t65 + (-t34 * t6 + t35 * t7) * mrSges(5,3) + (-Ifges(5,2) * t35 + t13 - t32) * t34 / 0.2e1;];
tau = t4;
