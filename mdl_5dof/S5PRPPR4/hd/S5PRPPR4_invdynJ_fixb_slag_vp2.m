% Calculate vector of inverse dynamics joint torques for
% S5PRPPR4
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
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-31 17:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRPPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:36:46
% EndTime: 2019-12-31 17:36:50
% DurationCPUTime: 2.85s
% Computational Cost: add. (850->210), mult. (1781->284), div. (0->0), fcn. (1118->6), ass. (0->93)
t74 = sin(pkin(8));
t75 = cos(pkin(8));
t124 = t74 ^ 2 + t75 ^ 2;
t56 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t113 = m(5) + m(6);
t94 = m(4) + t113;
t123 = Ifges(4,4) - Ifges(5,5);
t108 = t74 * mrSges(5,3);
t83 = t75 * mrSges(5,1) + t108;
t84 = mrSges(4,1) * t75 - mrSges(4,2) * t74;
t122 = -t84 - t83 - mrSges(3,1);
t121 = m(6) * pkin(6) + mrSges(3,2) + mrSges(6,3);
t103 = t74 * qJ(4);
t87 = pkin(2) + t103;
t42 = (pkin(3) + pkin(4)) * t75 + t87;
t20 = qJD(2) * t42 - qJD(3);
t76 = sin(qJ(5));
t100 = qJD(2) * t76;
t77 = cos(qJ(5));
t99 = qJD(2) * t77;
t40 = -t100 * t74 - t75 * t99;
t41 = -t100 * t75 + t74 * t99;
t120 = m(6) * t20 - mrSges(6,1) * t40 + mrSges(6,2) * t41 + t83 * qJD(2);
t73 = pkin(7) + qJ(2);
t69 = sin(t73);
t70 = cos(t73);
t119 = -g(1) * t70 - g(2) * t69;
t35 = t74 * qJDD(1) + t75 * t56;
t110 = t35 * t75;
t118 = t124 * t56 + t110 + t119;
t51 = -pkin(3) * t75 - t87;
t117 = t51 * qJDD(2);
t115 = t41 / 0.2e1;
t114 = m(2) + m(3);
t111 = Ifges(6,4) * t41;
t109 = t70 * t75;
t107 = -pkin(6) + qJ(3);
t97 = qJ(3) * qJD(2);
t48 = t74 * qJD(1) + t75 * t97;
t106 = t48 * t75 * qJD(3) + qJ(3) * t110;
t105 = t70 * pkin(2) + t69 * qJ(3);
t102 = qJD(2) * t74;
t101 = qJD(2) * t75;
t98 = qJD(4) * t74;
t96 = qJDD(2) * t74;
t95 = qJDD(2) * t75;
t86 = pkin(3) * t109 + t70 * t103 + t105;
t47 = qJD(1) * t75 - t74 * t97;
t34 = qJDD(1) * t75 - t74 * t56;
t80 = t74 * t76 + t75 * t77;
t38 = t80 * qJD(5);
t45 = t74 * t77 - t75 * t76;
t12 = -qJD(2) * t38 + qJDD(2) * t45;
t39 = t45 * qJD(5);
t13 = -qJD(2) * t39 - qJDD(2) * t80;
t82 = -t13 * mrSges(6,1) + t12 * mrSges(6,2);
t81 = -mrSges(6,1) * t80 - mrSges(6,2) * t45;
t43 = qJD(4) - t47;
t29 = -pkin(6) * t102 + t43;
t31 = -pkin(6) * t101 + t48;
t5 = t29 * t77 - t31 * t76;
t6 = t29 * t76 + t31 * t77;
t52 = t107 * t74;
t53 = t107 * t75;
t15 = t52 * t77 - t53 * t76;
t16 = t52 * t76 + t53 * t77;
t30 = qJDD(4) - t34;
t79 = -qJD(2) * t98 + qJDD(3);
t66 = -qJDD(2) * pkin(2) + qJDD(3);
t63 = mrSges(4,2) * t96;
t37 = t48 * t101;
t33 = qJD(2) * t51 + qJD(3);
t32 = Ifges(6,4) * t40;
t28 = t35 * t74;
t27 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t41;
t26 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t40;
t24 = t80 * t70;
t23 = t45 * t70;
t22 = t80 * t69;
t21 = t45 * t69;
t19 = -pkin(6) * t95 + t35;
t18 = t79 + t117;
t17 = -pkin(6) * t96 + t30;
t14 = qJDD(2) * t42 - t79;
t10 = Ifges(6,1) * t41 + Ifges(6,5) * qJD(5) + t32;
t9 = Ifges(6,2) * t40 + Ifges(6,6) * qJD(5) + t111;
t8 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t13;
t7 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t12;
t4 = qJD(3) * t45 - qJD(5) * t16;
t3 = qJD(3) * t80 + qJD(5) * t15;
t2 = -qJD(5) * t6 + t17 * t77 - t19 * t76;
t1 = qJD(5) * t5 + t17 * t76 + t19 * t77;
t11 = [-t38 * t26 - t39 * t27 - t80 * t7 + t45 * t8 + t114 * qJDD(1) + m(4) * (t34 * t75 + t28) + m(5) * (-t30 * t75 + t28) + m(6) * (t1 * t45 - t2 * t80 - t38 * t6 - t39 * t5) + (-t94 - t114) * g(3); m(5) * (t18 * t51 + (qJ(3) * t30 + qJD(3) * t43 - qJD(4) * t33) * t74 + t106) + m(4) * (-pkin(2) * t66 + (-qJ(3) * t34 - qJD(3) * t47) * t74 + t106) + (-mrSges(6,3) * t1 - Ifges(6,4) * t12 - Ifges(6,2) * t13 - Ifges(6,6) * qJDD(5)) * t80 - t66 * t84 - t14 * t81 + t42 * t82 + (-mrSges(6,3) * t2 + Ifges(6,1) * t12 + Ifges(6,4) * t13 + Ifges(6,5) * qJDD(5)) * t45 + (t22 * mrSges(6,1) + t21 * mrSges(6,2) + (m(4) * pkin(2) - m(5) * t51 + m(6) * t42 - t122) * t69 + (-t94 * qJ(3) + t121) * t70) * g(1) + m(6) * (t1 * t16 + t14 * t42 + t15 * t2 + t3 * t6 + t4 * t5) + (-m(4) * t105 - m(6) * (pkin(4) * t109 + t86) - t24 * mrSges(6,1) - t23 * mrSges(6,2) - m(5) * t86 + t121 * t69 + t122 * t70) * g(2) + (pkin(2) * mrSges(4,1) + (Ifges(5,3) + Ifges(4,2)) * t75 + t123 * t74) * t95 + (t123 * t75 + (Ifges(4,1) + Ifges(5,1)) * t74) * t96 + t120 * t98 + (-t34 * t74 + t118) * mrSges(4,3) + (t30 * t74 + t118) * mrSges(5,2) + (-t117 - t18) * t83 + t20 * (mrSges(6,1) * t39 - mrSges(6,2) * t38) + qJD(5) * (-Ifges(6,5) * t38 - Ifges(6,6) * t39) / 0.2e1 + t40 * (-Ifges(6,4) * t38 - Ifges(6,2) * t39) / 0.2e1 + (-Ifges(6,1) * t38 - Ifges(6,4) * t39) * t115 - pkin(2) * t63 - t38 * t10 / 0.2e1 - t39 * t9 / 0.2e1 + t3 * t26 + t4 * t27 + t15 * t7 + t16 * t8 + (t5 * t38 - t6 * t39) * mrSges(6,3) + Ifges(3,3) * qJDD(2); t40 * t26 - t41 * t27 + t63 + (-t108 + (-mrSges(4,1) - mrSges(5,1)) * t75) * qJDD(2) - t82 + (-g(1) * t69 + g(2) * t70) * t94 - (mrSges(5,2) + mrSges(4,3)) * qJD(2) ^ 2 * t124 + (t40 * t6 - t41 * t5 - t14) * m(6) + (-t102 * t43 + t18 - t37) * m(5) + (t102 * t47 - t37 + t66) * m(4); t77 * t7 + t76 * t8 + (t26 * t77 - t27 * t76) * qJD(5) + t113 * t75 * g(3) + m(5) * t30 + m(6) * (t1 * t76 + t2 * t77 + (-t5 * t76 + t6 * t77) * qJD(5)) + (qJDD(2) * mrSges(5,2) + (m(5) * t33 - t120) * qJD(2) + t113 * t119) * t74; Ifges(6,5) * t12 + Ifges(6,6) * t13 + Ifges(6,3) * qJDD(5) - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t20 * (mrSges(6,1) * t41 + mrSges(6,2) * t40) - t41 * (Ifges(6,1) * t40 - t111) / 0.2e1 + t9 * t115 - qJD(5) * (Ifges(6,5) * t40 - Ifges(6,6) * t41) / 0.2e1 - t5 * t26 + t6 * t27 - g(1) * (mrSges(6,1) * t23 - mrSges(6,2) * t24) - g(2) * (mrSges(6,1) * t21 - mrSges(6,2) * t22) - g(3) * t81 + (t40 * t5 + t41 * t6) * mrSges(6,3) - (-Ifges(6,2) * t41 + t10 + t32) * t40 / 0.2e1;];
tau = t11;
