% Calculate vector of inverse dynamics joint torques for
% S4RPRP6
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-12-31 16:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RPRP6_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP6_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:59
% EndTime: 2019-12-31 16:46:03
% DurationCPUTime: 2.01s
% Computational Cost: add. (527->180), mult. (1034->221), div. (0->0), fcn. (415->4), ass. (0->80)
t126 = Ifges(5,4) + Ifges(4,4);
t134 = Ifges(4,1) + Ifges(5,1);
t133 = Ifges(4,2) + Ifges(5,2);
t50 = cos(qJ(3));
t136 = t126 * t50;
t48 = sin(qJ(3));
t135 = t126 * t48;
t109 = t50 / 0.2e1;
t125 = Ifges(5,5) + Ifges(4,5);
t124 = Ifges(5,6) + Ifges(4,6);
t132 = -t133 * t48 + t136;
t131 = t134 * t50 - t135;
t40 = pkin(3) * t48 + qJ(2);
t28 = qJD(1) * t40 + qJD(4);
t130 = t28 * (mrSges(5,1) * t50 - mrSges(5,2) * t48) + (-t124 * t50 - t125 * t48) * qJD(3) / 0.2e1;
t129 = t48 / 0.2e1;
t52 = -pkin(1) - pkin(5);
t37 = qJD(1) * t52 + qJD(2);
t127 = -qJ(4) * qJD(1) + t37;
t51 = cos(qJ(1));
t107 = g(2) * t51;
t49 = sin(qJ(1));
t72 = -g(1) * t49 + t107;
t123 = t132 * qJD(1) + t124 * qJD(3);
t122 = t131 * qJD(1) + t125 * qJD(3);
t83 = qJD(1) * qJD(2);
t38 = qJDD(1) * qJ(2) + t83;
t90 = qJD(1) * t48;
t32 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t90;
t89 = qJD(1) * t50;
t34 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t89;
t120 = t50 * t32 - t48 * t34;
t36 = qJDD(1) * t52 + qJDD(2);
t88 = qJD(3) * t48;
t3 = t50 * t36 - t37 * t88;
t87 = qJD(3) * t50;
t4 = t48 * t36 + t37 * t87;
t71 = t3 * t50 + t4 * t48;
t82 = qJD(1) * qJD(3);
t25 = qJDD(1) * t50 - t48 * t82;
t81 = qJD(1) * qJD(4);
t1 = qJDD(3) * pkin(3) - qJ(4) * t25 - t50 * t81 + t3;
t26 = -qJDD(1) * t48 - t50 * t82;
t2 = qJ(4) * t26 - t48 * t81 + t4;
t119 = t1 * t50 + t2 * t48;
t118 = -m(5) - m(4) - m(3);
t117 = mrSges(2,1) + mrSges(4,3) - mrSges(3,2);
t69 = mrSges(4,1) * t48 + mrSges(4,2) * t50;
t78 = m(5) * pkin(3) + mrSges(5,1);
t97 = t50 * mrSges(5,2);
t116 = -t48 * t78 + mrSges(2,2) - mrSges(3,3) - t69 - t97;
t114 = -qJ(2) * (mrSges(4,1) * t50 - mrSges(4,2) * t48) + (-t133 * t50 - t135) * t129 - (-t134 * t48 - t136) * t50 / 0.2e1;
t91 = qJ(4) - t52;
t84 = qJDD(1) * mrSges(3,2);
t75 = -t26 * mrSges(5,1) + t25 * mrSges(5,2);
t30 = t91 * t50;
t74 = (t38 + t83) * qJ(2);
t10 = t127 * t50;
t6 = qJD(3) * pkin(3) + t10;
t9 = t127 * t48;
t70 = -t6 * t48 + t9 * t50;
t68 = mrSges(5,1) * t48 + t97;
t53 = qJD(1) ^ 2;
t55 = -t53 * qJ(2) + t72;
t47 = -qJ(4) - pkin(5);
t42 = -qJDD(1) * pkin(1) + qJDD(2);
t39 = pkin(3) * t87 + qJD(2);
t33 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t89;
t31 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t90;
t29 = t91 * t48;
t22 = t69 * qJD(1);
t21 = t68 * qJD(1);
t14 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t26;
t13 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t26;
t12 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t25;
t11 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t25;
t8 = -qJD(3) * t30 - qJD(4) * t48;
t7 = -qJD(4) * t50 + t88 * t91;
t5 = -pkin(3) * t26 + qJDD(4) + t38;
t15 = [((-m(5) * (-pkin(1) + t47) + mrSges(5,3) + m(3) * pkin(1) - m(4) * t52 + t117) * t49 + (t118 * qJ(2) + t116) * t51) * g(1) + (t6 * t88 - t87 * t9 - t107 - t119) * mrSges(5,3) - t123 * t87 / 0.2e1 - t71 * mrSges(4,3) - t122 * t88 / 0.2e1 + (t118 * (t51 * pkin(1) + t49 * qJ(2)) + (-m(4) * pkin(5) + m(5) * t47 - t117) * t51 + t116 * t49) * g(2) - t114 * t82 + (t120 * t52 + t130) * qJD(3) + (Ifges(3,1) + Ifges(2,3)) * qJDD(1) + (0.2e1 * t125 * t109 - t124 * t48) * qJDD(3) + t40 * t75 + m(3) * (-pkin(1) * t42 + t74) + (t71 * m(4) + t50 * t12 + t48 * t14) * t52 + t5 * t68 + m(5) * (-t1 * t30 - t2 * t29 + t28 * t39 + t40 * t5 + t6 * t7 + t8 * t9) + m(4) * t74 - (t126 * t25 + t133 * t26) * t48 / 0.2e1 + (t126 * t26 + t134 * t25) * t109 + t42 * mrSges(3,2) + qJ(2) * (-t26 * mrSges(4,1) + t25 * mrSges(4,2)) - t29 * t13 - t30 * t11 + t8 * t31 + t7 * t33 + t39 * t21 + qJD(2) * t22 + t131 * t25 / 0.2e1 + t132 * t26 / 0.2e1 - pkin(1) * t84 + (t69 + 0.2e1 * mrSges(3,3)) * t38; t84 - t53 * mrSges(3,3) + (t11 + t12) * t50 + (t13 + t14) * t48 + (-t21 - t22) * qJD(1) + ((t31 + t32) * t50 + (-t33 - t34) * t48) * qJD(3) + (-t28 * qJD(1) + qJD(3) * t70 + t119 + t72) * m(5) + (t55 + t71) * m(4) + (t42 + t55) * m(3); t3 * mrSges(4,1) + t1 * mrSges(5,1) - t4 * mrSges(4,2) - t2 * mrSges(5,2) - t10 * t31 - t120 * t37 + t124 * t26 + t125 * t25 + (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + (t69 + t68) * g(3) + (t11 + (g(3) * t48 + t1) * m(5)) * pkin(3) + (t33 - m(5) * (t10 - t6)) * t9 + (t114 * qJD(1) + (-m(5) * t28 - t21) * t50 * pkin(3) + t70 * mrSges(5,3) + t122 * t129 + t123 * t109 - t130) * qJD(1) - t72 * ((mrSges(4,2) + mrSges(5,2)) * t48 - t50 * (mrSges(4,1) + t78)); (t48 * t31 + t50 * t33) * qJD(1) + (t5 - g(1) * t51 - g(2) * t49 - (-t48 * t9 - t50 * t6) * qJD(1)) * m(5) + t75;];
tau = t15;
