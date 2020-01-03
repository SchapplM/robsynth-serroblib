% Calculate vector of inverse dynamics joint torques for
% S4RPRP4
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-12-31 16:44
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:43:38
% EndTime: 2019-12-31 16:43:42
% DurationCPUTime: 2.19s
% Computational Cost: add. (589->205), mult. (1232->264), div. (0->0), fcn. (570->8), ass. (0->84)
t126 = Ifges(5,4) + Ifges(4,5);
t125 = Ifges(5,6) - Ifges(4,6);
t51 = cos(pkin(6));
t111 = pkin(1) * t51;
t40 = -pkin(2) - t111;
t52 = sin(qJ(3));
t54 = cos(qJ(3));
t71 = t52 * mrSges(5,1) - t54 * mrSges(5,3);
t73 = mrSges(4,1) * t52 + mrSges(4,2) * t54;
t66 = pkin(3) * t54 + qJ(4) * t52;
t63 = -pkin(2) - t66;
t20 = t63 - t111;
t8 = t20 * qJD(1);
t130 = -t40 * qJD(1) * t73 - t8 * t71;
t49 = qJ(1) + pkin(6);
t44 = sin(t49);
t45 = cos(t49);
t129 = g(1) * t45 + g(2) * t44;
t128 = -m(5) - m(4);
t74 = mrSges(4,1) * t54 - mrSges(4,2) * t52;
t127 = -mrSges(3,1) - t74;
t90 = qJD(1) * qJD(3);
t25 = -qJDD(1) * t54 + t52 * t90;
t14 = -mrSges(5,2) * t25 + qJDD(3) * mrSges(5,3);
t124 = -qJDD(3) * mrSges(4,2) - mrSges(4,3) * t25 + t14;
t26 = qJDD(1) * t52 + t54 * t90;
t13 = -qJDD(3) * mrSges(5,1) + mrSges(5,2) * t26;
t123 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t26 - t13;
t91 = t54 * qJD(1);
t43 = Ifges(4,4) * t91;
t100 = Ifges(5,5) * t54;
t70 = Ifges(5,1) * t52 - t100;
t92 = t52 * qJD(1);
t122 = Ifges(4,1) * t92 + qJD(1) * t70 + qJD(3) * t126 + t43;
t82 = mrSges(4,3) * t92;
t83 = mrSges(5,2) * t92;
t121 = t82 + t83 + (-mrSges(4,1) - mrSges(5,1)) * qJD(3);
t81 = mrSges(5,2) * t91;
t34 = qJD(3) * mrSges(5,3) + t81;
t86 = mrSges(4,3) * t91;
t120 = -qJD(3) * mrSges(4,2) + t34 + t86;
t118 = t125 * t52 + t126 * t54;
t50 = sin(pkin(6));
t39 = pkin(1) * t50 + pkin(5);
t29 = t39 * qJD(1);
t27 = t39 * qJDD(1);
t93 = qJD(3) * t54;
t88 = qJD(2) * t93 + qJDD(2) * t52 + t27 * t54;
t94 = qJD(3) * t52;
t4 = -t29 * t94 + t88;
t10 = qJD(2) * t52 + t29 * t54;
t5 = -qJD(3) * t10 + qJDD(2) * t54 - t27 * t52;
t117 = t4 * t54 - t5 * t52;
t99 = t29 * t52;
t1 = qJDD(3) * qJ(4) + (qJD(4) - t99) * qJD(3) + t88;
t2 = -qJDD(3) * pkin(3) + qJDD(4) - t5;
t116 = t1 * t54 + t2 * t52;
t112 = t52 / 0.2e1;
t53 = sin(qJ(1));
t110 = pkin(1) * t53;
t55 = cos(qJ(1));
t48 = t55 * pkin(1);
t103 = Ifges(4,4) * t52;
t102 = Ifges(4,4) * t54;
t101 = Ifges(5,5) * t52;
t76 = -t90 / 0.2e1;
t72 = t54 * mrSges(5,1) + t52 * mrSges(5,3);
t69 = Ifges(4,2) * t54 + t103;
t65 = pkin(3) * t52 - qJ(4) * t54;
t28 = t40 * qJDD(1);
t9 = qJD(2) * t54 - t99;
t60 = t52 * (Ifges(4,1) * t54 - t103);
t59 = t54 * (Ifges(5,3) * t52 + t100);
t58 = m(5) * t66 + t72;
t42 = Ifges(5,5) * t92;
t24 = t65 * qJD(1);
t23 = t72 * qJD(1);
t17 = Ifges(4,6) * qJD(3) + qJD(1) * t69;
t16 = Ifges(5,6) * qJD(3) - Ifges(5,3) * t91 + t42;
t15 = qJD(3) * t65 - qJD(4) * t52;
t7 = qJD(3) * qJ(4) + t10;
t6 = -qJD(3) * pkin(3) + qJD(4) - t9;
t3 = pkin(3) * t25 - qJ(4) * t26 - qJD(4) * t92 + t28;
t11 = [-t54 * (Ifges(5,5) * t26 + Ifges(5,6) * qJDD(3)) / 0.2e1 + t54 * (Ifges(4,4) * t26 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t122 * t93 / 0.2e1 + (-t125 * t54 + t126 * t52) * qJDD(3) / 0.2e1 + ((Ifges(5,1) + Ifges(4,1)) * t26 + t126 * qJDD(3)) * t112 + (t6 * t93 - t7 * t94 + t116 - t129) * mrSges(5,2) + (-t10 * t94 - t9 * t93 + t117 - t129) * mrSges(4,3) - t20 * t26 * mrSges(5,3) + t40 * t26 * mrSges(4,2) + (-t69 / 0.2e1 + t40 * mrSges(4,1) + t20 * mrSges(5,1) + t101 / 0.2e1 + (-Ifges(5,3) - Ifges(4,2) / 0.2e1) * t54 + (Ifges(5,5) - Ifges(4,4)) * t112) * t25 + m(5) * (t15 * t8 + t20 * t3) + (-t17 / 0.2e1 + t16 / 0.2e1) * t94 - t15 * t23 - t3 * t72 + (t52 * (Ifges(5,1) * t54 + t101) + t54 * (-Ifges(4,2) * t52 + t102) + t60) * t90 / 0.2e1 + (Ifges(4,1) * t52 + t102 + t70) * t26 / 0.2e1 + (m(3) * t110 + mrSges(2,1) * t53 + mrSges(2,2) * t55 + mrSges(3,2) * t45 + t128 * (pkin(5) * t45 - t110) + (m(4) * pkin(2) - m(5) * t63 - t127 + t72) * t44) * g(1) + (-m(3) * t48 - mrSges(2,1) * t55 + mrSges(2,2) * t53 + mrSges(3,2) * t44 + t128 * (pkin(2) * t45 + pkin(5) * t44 + t48) + (-t58 + t127) * t45) * g(2) + t59 * t76 + (m(4) * t40 - t74) * t28 + (Ifges(3,3) + Ifges(2,3) + (0.2e1 * t51 * mrSges(3,1) - 0.2e1 * t50 * mrSges(3,2) + m(3) * (t50 ^ 2 + t51 ^ 2) * pkin(1)) * pkin(1)) * qJDD(1) + (t118 * qJD(3) / 0.2e1 - t130) * qJD(3) + (m(5) * ((-t52 * t7 + t54 * t6) * qJD(3) + t116) + m(4) * ((-t10 * t52 - t54 * t9) * qJD(3) + t117) - t120 * t94 + t121 * t93 - t123 * t52 + t124 * t54) * t39; m(3) * qJDD(2) + t123 * t54 + t124 * t52 + (t120 * t54 + t121 * t52) * qJD(3) + m(4) * (t4 * t52 + t5 * t54 + (t10 * t54 - t52 * t9) * qJD(3)) + m(5) * (t1 * t52 - t2 * t54 + (t52 * t6 + t54 * t7) * qJD(3)) + (-m(3) + t128) * g(3); (-m(5) * t7 - t120 + t86) * t9 + (-m(5) * t6 - t121 + t82) * t10 - (-Ifges(4,2) * t92 + t122 + t43) * t91 / 0.2e1 + t125 * t25 + t126 * t26 + t118 * t76 + t7 * t83 + (-pkin(3) * t2 + qJ(4) * t1 + qJD(4) * t7 - t24 * t8) * m(5) + (Ifges(5,2) + Ifges(4,3)) * qJDD(3) + qJD(4) * t34 - pkin(3) * t13 + qJ(4) * t14 + t24 * t23 + t1 * mrSges(5,3) - t2 * mrSges(5,1) - t4 * mrSges(4,2) + t5 * mrSges(4,1) - (Ifges(5,1) * t91 + t16 + t42) * t92 / 0.2e1 + (-t58 - t74) * g(3) + t17 * t92 / 0.2e1 - t6 * t81 + t129 * (m(5) * t65 + t71 + t73) + ((-t60 / 0.2e1 + t59 / 0.2e1) * qJD(1) + t130) * qJD(1); -t23 * t92 - qJD(3) * t34 + (g(3) * t54 - t7 * qJD(3) - t129 * t52 + t8 * t92 + t2) * m(5) + t13;];
tau = t11;
