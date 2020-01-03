% Calculate vector of inverse dynamics joint torques for
% S4PRPR4
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
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-12-31 16:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:21:56
% EndTime: 2019-12-31 16:21:58
% DurationCPUTime: 0.72s
% Computational Cost: add. (336->114), mult. (575->158), div. (0->0), fcn. (230->4), ass. (0->53)
t26 = sin(qJ(4));
t71 = -t26 / 0.2e1;
t27 = cos(qJ(4));
t62 = t27 / 0.2e1;
t25 = pkin(6) + qJ(2);
t22 = sin(t25);
t23 = cos(t25);
t65 = -g(1) * t22 + g(2) * t23;
t70 = -m(5) - m(4);
t69 = mrSges(4,2) - mrSges(3,1);
t38 = mrSges(5,1) * t26 + mrSges(5,2) * t27;
t68 = -t38 + mrSges(3,2);
t44 = qJD(2) * qJD(3);
t17 = qJDD(2) * qJ(3) + t44;
t43 = qJD(2) * qJD(4);
t10 = qJDD(2) * t27 - t26 * t43;
t3 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t10;
t11 = -qJDD(2) * t26 - t27 * t43;
t4 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t11;
t67 = t26 * t4 + t27 * t3;
t28 = -pkin(2) - pkin(5);
t15 = t28 * qJDD(2) + qJDD(3);
t16 = t28 * qJD(2) + qJD(3);
t5 = -qJD(1) * t26 + t16 * t27;
t1 = t5 * qJD(4) + qJDD(1) * t27 + t15 * t26;
t6 = qJD(1) * t27 + t16 * t26;
t2 = -t6 * qJD(4) - qJDD(1) * t26 + t15 * t27;
t66 = t1 * t26 + t2 * t27;
t50 = mrSges(5,3) * qJD(2);
t12 = -qJD(4) * mrSges(5,2) - t26 * t50;
t13 = qJD(4) * mrSges(5,1) - t27 * t50;
t64 = (t27 * t12 - t26 * t13) * qJD(4);
t39 = mrSges(5,1) * t27 - mrSges(5,2) * t26;
t54 = Ifges(5,4) * t27;
t55 = Ifges(5,4) * t26;
t63 = (-Ifges(5,2) * t27 - t55) * t71 + (-Ifges(5,1) * t26 - t54) * t62 + qJ(3) * t39;
t49 = qJD(4) * t26;
t48 = qJD(4) * t27;
t46 = qJDD(2) * mrSges(4,2);
t45 = m(2) + m(3) + m(4);
t41 = (t17 + t44) * qJ(3);
t40 = -t26 * t5 + t27 * t6;
t37 = Ifges(5,1) * t27 - t55;
t36 = -Ifges(5,2) * t26 + t54;
t35 = -Ifges(5,5) * t26 - Ifges(5,6) * t27;
t29 = qJD(2) ^ 2;
t31 = -t29 * qJ(3) + t65;
t30 = t40 * qJD(4) + t66;
t21 = -qJDD(2) * pkin(2) + qJDD(3);
t9 = t38 * qJD(2);
t8 = Ifges(5,5) * qJD(4) + t37 * qJD(2);
t7 = Ifges(5,6) * qJD(4) + t36 * qJD(2);
t14 = [-t12 * t49 - t13 * t48 + t27 * t4 - t26 * t3 + m(5) * (t1 * t27 - t2 * t26 + (-t26 * t6 - t27 * t5) * qJD(4)) + t45 * qJDD(1) + (-m(5) - t45) * g(3); -pkin(2) * t46 - t7 * t48 / 0.2e1 - t8 * t49 / 0.2e1 + m(5) * t41 + m(4) * (-pkin(2) * t21 + t41) + qJD(4) ^ 2 * t35 / 0.2e1 + t11 * t36 / 0.2e1 + t10 * t37 / 0.2e1 + (Ifges(5,4) * t10 + Ifges(5,2) * t11) * t71 + qJ(3) * (-t11 * mrSges(5,1) + t10 * mrSges(5,2)) + t17 * t38 + t21 * mrSges(4,2) + qJD(3) * t9 + (Ifges(5,1) * t10 + Ifges(5,4) * t11) * t62 + t28 * t64 + (t30 * m(5) + t67) * t28 + (Ifges(4,1) + Ifges(3,3)) * qJDD(2) + (0.2e1 * Ifges(5,5) * t62 - Ifges(5,6) * t26) * qJDD(4) + t63 * t43 + (t70 * (t23 * pkin(2) + t22 * qJ(3)) + (-m(5) * pkin(5) + t69) * t23 + t68 * t22) * g(2) + ((m(4) * pkin(2) - m(5) * t28 - t69) * t22 + (t70 * qJ(3) + t68) * t23) * g(1) + (-t6 * t48 + t5 * t49 - t65 - t66) * mrSges(5,3) + (-g(1) * t23 - g(2) * t22 + 0.2e1 * t17) * mrSges(4,3); t46 - t29 * mrSges(4,3) - qJD(2) * t9 + t64 + (t30 + t31) * m(5) + (t21 + t31) * m(4) + t67; t2 * mrSges(5,1) - t1 * mrSges(5,2) + Ifges(5,5) * t10 + Ifges(5,6) * t11 + Ifges(5,3) * qJDD(4) + g(3) * t38 - t5 * t12 + t6 * t13 + (t26 * t8 / 0.2e1 + t7 * t62 - qJD(4) * t35 / 0.2e1 - t63 * qJD(2) + t40 * mrSges(5,3)) * qJD(2) + t65 * t39;];
tau = t14;
