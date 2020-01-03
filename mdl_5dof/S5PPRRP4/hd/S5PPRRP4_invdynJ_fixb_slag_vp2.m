% Calculate vector of inverse dynamics joint torques for
% S5PPRRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:31
% DurationCPUTime: 2.41s
% Computational Cost: add. (745->244), mult. (1650->305), div. (0->0), fcn. (901->6), ass. (0->100)
t134 = Ifges(6,4) + Ifges(5,4);
t58 = sin(qJ(4));
t60 = cos(qJ(4));
t128 = -mrSges(6,1) * t60 + t58 * mrSges(6,2);
t41 = -t60 * mrSges(5,1) + t58 * mrSges(5,2);
t139 = -t41 - t128;
t101 = qJD(1) * t60;
t59 = sin(qJ(3));
t100 = qJD(2) * t59;
t43 = qJD(3) * pkin(6) + t100;
t16 = -t43 * t58 - t101;
t55 = t58 * qJD(1);
t17 = t43 * t60 - t55;
t61 = cos(qJ(3));
t93 = qJD(2) * qJD(3);
t50 = t61 * t93;
t35 = t59 * qJDD(2) + t50;
t26 = qJDD(3) * pkin(6) + t35;
t3 = t16 * qJD(4) - qJDD(1) * t58 + t60 * t26;
t65 = qJD(4) * t55 + (-qJD(4) * t43 - qJDD(1)) * t60;
t4 = -t26 * t58 + t65;
t72 = t3 * t60 - t4 * t58;
t97 = qJD(4) * t60;
t98 = qJD(4) * t58;
t138 = -t16 * t97 - t17 * t98 + t72;
t137 = m(3) + m(4);
t136 = -m(5) - m(6);
t135 = Ifges(5,1) + Ifges(6,1);
t133 = Ifges(6,5) + Ifges(5,5);
t132 = Ifges(6,2) + Ifges(5,2);
t131 = Ifges(6,6) + Ifges(5,6);
t30 = t128 * qJD(3);
t130 = -t41 * qJD(3) - t30;
t127 = t134 * t60;
t51 = pkin(4) * t60 + pkin(3);
t126 = m(5) * pkin(3) + m(6) * t51 + mrSges(4,1) + t139;
t57 = -qJ(5) - pkin(6);
t125 = m(5) * pkin(6) - m(6) * t57 - mrSges(4,2) + mrSges(5,3) + mrSges(6,3);
t92 = qJD(3) * qJD(4);
t33 = qJDD(3) * t58 + t60 * t92;
t91 = qJD(3) * qJD(5);
t1 = qJDD(4) * pkin(4) - qJ(5) * t33 + (-t26 - t91) * t58 + t65;
t75 = qJ(5) * qJD(3) + t43;
t10 = t60 * t75 - t55;
t32 = qJDD(3) * t60 - t58 * t92;
t2 = qJ(5) * t32 + t60 * t91 + t3;
t9 = -t58 * t75 - t101;
t8 = qJD(4) * pkin(4) + t9;
t123 = -t1 * t58 - t10 * t98 + t2 * t60 - t8 * t97;
t95 = t60 * qJD(3);
t38 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t95;
t39 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t95;
t105 = t38 + t39;
t96 = t58 * qJD(3);
t36 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t96;
t37 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t96;
t106 = t36 + t37;
t122 = t105 * t60 - t106 * t58;
t112 = Ifges(5,4) * t58;
t110 = Ifges(6,4) * t58;
t12 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t32;
t13 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t32;
t108 = t12 + t13;
t14 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t33;
t15 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t33;
t107 = -t14 - t15;
t103 = cos(pkin(7));
t102 = sin(pkin(7));
t99 = qJD(2) * t61;
t94 = m(2) + t137;
t89 = pkin(4) * t98;
t87 = t58 * t99;
t86 = t60 * t99;
t49 = t59 * t93;
t6 = -t32 * mrSges(6,1) + t33 * mrSges(6,2);
t77 = qJD(4) * t57;
t34 = qJDD(2) * t61 - t49;
t71 = t10 * t60 - t58 * t8;
t69 = Ifges(5,2) * t60 + t112;
t68 = Ifges(6,2) * t60 + t110;
t67 = t16 * t58 - t17 * t60;
t25 = -qJDD(3) * pkin(3) - t34;
t62 = qJD(3) ^ 2;
t53 = Ifges(5,4) * t95;
t52 = Ifges(6,4) * t95;
t44 = -qJD(3) * pkin(3) - t99;
t42 = t57 * t60;
t40 = t57 * t58;
t29 = -t102 * t61 + t103 * t59;
t28 = -t102 * t59 - t103 * t61;
t24 = Ifges(5,1) * t96 + Ifges(5,5) * qJD(4) + t53;
t23 = Ifges(6,1) * t96 + Ifges(6,5) * qJD(4) + t52;
t22 = Ifges(5,6) * qJD(4) + qJD(3) * t69;
t21 = Ifges(6,6) * qJD(4) + qJD(3) * t68;
t20 = -qJD(3) * t51 + qJD(5) - t99;
t19 = -qJD(5) * t58 + t60 * t77;
t18 = qJD(5) * t60 + t58 * t77;
t7 = -mrSges(5,1) * t32 + mrSges(5,2) * t33;
t5 = -pkin(4) * t32 + qJDD(5) + t25;
t11 = [t107 * t60 - t108 * t58 - t122 * qJD(4) + m(5) * (qJD(4) * t67 - t3 * t58 - t4 * t60) + m(6) * (-qJD(4) * t71 - t1 * t60 - t2 * t58) + t94 * qJDD(1) + (-t94 + t136) * g(3); m(3) * qJDD(2) + (qJDD(3) * mrSges(4,1) - t62 * mrSges(4,2) - t6 - t7 + t122 * qJD(3) + m(4) * t34 + m(6) * (t10 * t95 - t8 * t96 - t5) + m(5) * (-t16 * t96 + t17 * t95 - t25)) * t61 + (-t62 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t108 * t60 + t107 * t58 - t130 * qJD(3) + (-t105 * t58 - t106 * t60) * qJD(4) + m(4) * t35 + m(6) * (qJD(3) * t20 + t123) + m(5) * (qJD(3) * t44 + t138)) * t59 + (-g(1) * t102 + g(2) * t103) * (-t136 + t137); (t19 + t87) * t36 + (t20 * (mrSges(6,1) * t58 + mrSges(6,2) * t60) + t44 * (mrSges(5,1) * t58 + mrSges(5,2) * t60)) * qJD(4) + t127 * t33 / 0.2e1 - t39 * t86 + ((-t132 * t58 + t127) * t60 + (t135 * t60 - t110 - t112) * t58) * t92 / 0.2e1 + t5 * t128 + t130 * t100 + (t131 * t60 + t133 * t58) * qJDD(4) + (-t131 * t58 + t133 * t60) * qJD(4) ^ 2 / 0.2e1 + (t125 * t29 - t126 * t28) * g(2) + (t125 * t28 + t126 * t29) * g(1) + t123 * mrSges(6,3) + (t132 * t32 + t134 * t33) * t60 / 0.2e1 + (t34 + t49) * mrSges(4,1) + (-t35 + t50) * mrSges(4,2) + (-t86 + t18) * t38 + (-pkin(3) * t25 - (t44 * t59 - t61 * t67) * qJD(2)) * m(5) + t40 * t14 + t25 * t41 - t42 * t12 - t51 * t6 - pkin(3) * t7 + (t134 * t58 + t68 + t69) * t32 / 0.2e1 + t138 * mrSges(5,3) + (t24 + t23) * t97 / 0.2e1 - (t22 + t21) * t98 / 0.2e1 + (t1 * t40 + t10 * t18 + t19 * t8 - t2 * t42 + t20 * t89 - t5 * t51 - (t20 * t59 + t61 * t71) * qJD(2)) * m(6) + t37 * t87 + t30 * t89 + t135 * t58 * t33 + (-t37 * t97 - t39 * t98 + m(5) * ((-t16 * t60 - t17 * t58) * qJD(4) + t72) - t58 * t15 + t60 * t13) * pkin(6) + Ifges(4,3) * qJDD(3); t4 * mrSges(5,1) + t1 * mrSges(6,1) - t3 * mrSges(5,2) - t2 * mrSges(6,2) - t16 * t39 + t17 * t37 - t9 * t38 + t133 * t33 + t131 * t32 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + t139 * g(3) + (t14 + (g(3) * t60 + t1) * m(6)) * pkin(4) + (t36 - m(6) * (-t8 + t9)) * t10 + ((-t44 * mrSges(5,2) - t20 * mrSges(6,2) + t16 * mrSges(5,3) + t8 * mrSges(6,3) - t23 / 0.2e1 - t24 / 0.2e1 - t52 / 0.2e1 - t53 / 0.2e1 + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4)) * t60 + (-t44 * mrSges(5,1) - t20 * mrSges(6,1) + t17 * mrSges(5,3) + t10 * mrSges(6,3) + t21 / 0.2e1 + t22 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t96 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t20 - t30) * pkin(4) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1 - Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1) * t95) * t58) * qJD(3) + (g(1) * t28 + g(2) * t29) * ((-mrSges(5,2) - mrSges(6,2)) * t60 + (-m(6) * pkin(4) - mrSges(5,1) - mrSges(6,1)) * t58); (t58 * t36 - t60 * t38) * qJD(3) + (-g(1) * t29 + g(2) * t28 - qJD(3) * t71 + t5) * m(6) + t6;];
tau = t11;
