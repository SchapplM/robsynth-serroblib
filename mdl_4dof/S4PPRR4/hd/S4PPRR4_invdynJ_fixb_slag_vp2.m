% Calculate vector of inverse dynamics joint torques for
% S4PPRR4
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
%   pkin=[a2,a3,a4,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-31 16:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PPRR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:18:26
% EndTime: 2019-12-31 16:18:29
% DurationCPUTime: 0.94s
% Computational Cost: add. (441->134), mult. (1062->200), div. (0->0), fcn. (709->10), ass. (0->70)
t39 = sin(qJ(4));
t41 = cos(qJ(4));
t53 = -t41 * mrSges(5,1) + t39 * mrSges(5,2);
t97 = m(5) * pkin(3) + mrSges(4,1) - t53;
t96 = -m(5) * pkin(5) + mrSges(4,2) - mrSges(5,3);
t35 = sin(pkin(7));
t37 = cos(pkin(7));
t40 = sin(qJ(3));
t42 = cos(qJ(3));
t95 = -t35 * t40 + t42 * t37;
t21 = t35 * t42 + t40 * t37;
t87 = t39 / 0.2e1;
t94 = m(4) + m(3);
t65 = qJD(3) * t39;
t25 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t65;
t64 = qJD(3) * t41;
t26 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t64;
t49 = -t39 * t25 + t41 * t26;
t93 = -qJD(3) * mrSges(4,2) + t49;
t60 = qJD(3) * qJD(4);
t23 = qJDD(3) * t41 - t39 * t60;
t12 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t23;
t24 = qJDD(3) * t39 + t41 * t60;
t13 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t24;
t92 = t41 * t12 - t39 * t13;
t61 = qJD(1) * qJD(3);
t5 = t21 * qJDD(1) + t95 * t61;
t3 = qJDD(3) * pkin(5) + t5;
t15 = t21 * qJD(1);
t11 = qJD(3) * pkin(5) + t15;
t7 = qJD(2) * t41 - t11 * t39;
t1 = t7 * qJD(4) + qJDD(2) * t39 + t3 * t41;
t8 = qJD(2) * t39 + t11 * t41;
t2 = -t8 * qJD(4) + qJDD(2) * t41 - t3 * t39;
t91 = t1 * t41 - t2 * t39;
t14 = t95 * qJD(1);
t89 = t95 * qJDD(1);
t10 = -qJD(3) * pkin(3) - t14;
t52 = mrSges(5,1) * t39 + mrSges(5,2) * t41;
t88 = t10 * t52 + qJD(4) * (Ifges(5,5) * t41 - Ifges(5,6) * t39) / 0.2e1;
t81 = Ifges(5,4) * t39;
t80 = Ifges(5,4) * t41;
t79 = Ifges(5,2) * t41;
t36 = sin(pkin(6));
t76 = t36 * t39;
t75 = t36 * t41;
t38 = cos(pkin(6));
t74 = t38 * t39;
t73 = t38 * t41;
t63 = qJD(4) * t39;
t62 = qJD(4) * t41;
t57 = -g(1) * t36 + g(2) * t38;
t56 = t8 * t39 + t7 * t41;
t55 = -t39 * t7 + t41 * t8;
t51 = t79 + t81;
t46 = t39 * (Ifges(5,1) * t41 - t81);
t17 = t21 * qJD(3);
t44 = -t56 * qJD(4) + t91;
t34 = pkin(7) + qJ(3);
t33 = cos(t34);
t32 = sin(t34);
t31 = Ifges(5,4) * t64;
t22 = t53 * qJD(3);
t19 = Ifges(5,1) * t65 + Ifges(5,5) * qJD(4) + t31;
t18 = Ifges(5,6) * qJD(4) + t51 * qJD(3);
t16 = t95 * qJD(3);
t9 = -mrSges(5,1) * t23 + mrSges(5,2) * t24;
t6 = -t21 * t61 + t89;
t4 = -qJDD(3) * pkin(3) + qJD(1) * t17 - t89;
t20 = [t17 * t22 - t95 * t9 + (-qJD(3) * t17 + qJDD(3) * t95) * mrSges(4,1) + t93 * t16 + (-qJDD(3) * mrSges(4,2) + (-t41 * t25 - t39 * t26) * qJD(4) + t92) * t21 + m(4) * (-t14 * t17 + t15 * t16 + t21 * t5 + t6 * t95) + m(5) * (t10 * t17 + t55 * t16 + t44 * t21 - t4 * t95) + (-m(2) - m(5) - t94) * g(3) + (m(2) + m(3) * (t35 ^ 2 + t37 ^ 2)) * qJDD(1); t39 * t12 + t41 * t13 + t49 * qJD(4) + (t55 * qJD(4) + t1 * t39 + t2 * t41 + t57) * m(5) + t94 * (qJDD(2) + t57); t24 * (Ifges(5,1) * t39 + t80) / 0.2e1 + t19 * t62 / 0.2e1 - t18 * t63 / 0.2e1 + t4 * t53 + t23 * t51 / 0.2e1 + (Ifges(5,1) * t24 + Ifges(5,4) * t23) * t87 + Ifges(4,3) * qJDD(3) - t5 * mrSges(4,2) + t6 * mrSges(4,1) + t41 * (Ifges(5,4) * t24 + Ifges(5,2) * t23) / 0.2e1 + (t46 + t41 * (-Ifges(5,2) * t39 + t80)) * t60 / 0.2e1 + t88 * qJD(4) + (t96 * t32 - t97 * t33) * g(3) + (-m(5) * t4 - t9) * pkin(3) + (0.2e1 * Ifges(5,5) * t87 + Ifges(5,6) * t41) * qJDD(4) + (-m(5) * t10 + qJD(3) * mrSges(4,1) - t22) * t15 + (-t55 * m(5) - t93) * t14 + (m(5) * t44 - t25 * t62 - t26 * t63 + t92) * pkin(5) + (-t7 * t62 - t8 * t63 + t91) * mrSges(5,3) + (t97 * t32 + t96 * t33) * (g(1) * t38 + g(2) * t36); Ifges(5,5) * t24 + Ifges(5,6) * t23 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t7 * t26 + t8 * t25 - g(1) * ((-t33 * t74 + t75) * mrSges(5,1) + (-t33 * t73 - t76) * mrSges(5,2)) - g(2) * ((-t33 * t76 - t73) * mrSges(5,1) + (-t33 * t75 + t74) * mrSges(5,2)) + g(3) * t52 * t32 + (t18 * t87 + (-t46 / 0.2e1 + t79 * t87) * qJD(3) + t56 * mrSges(5,3) - (t19 + t31) * t41 / 0.2e1 - t88) * qJD(3);];
tau = t20;
