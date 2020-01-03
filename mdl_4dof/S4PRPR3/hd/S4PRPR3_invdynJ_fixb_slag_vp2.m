% Calculate vector of inverse dynamics joint torques for
% S4PRPR3
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
% Datum: 2019-12-31 16:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4PRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:20:48
% EndTime: 2019-12-31 16:20:51
% DurationCPUTime: 1.12s
% Computational Cost: add. (573->149), mult. (1284->209), div. (0->0), fcn. (829->8), ass. (0->67)
t38 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t55 = sin(pkin(7));
t56 = cos(pkin(7));
t90 = t55 ^ 2 + t56 ^ 2;
t44 = t56 * qJDD(1);
t23 = -t38 * t55 + t44;
t24 = t55 * qJDD(1) + t56 * t38;
t88 = -t23 * t55 + t24 * t56;
t41 = pkin(3) * t56 + pkin(2);
t53 = pkin(7) + qJ(4);
t47 = sin(t53);
t49 = cos(t53);
t65 = t49 * mrSges(5,1) - t47 * mrSges(5,2);
t66 = -t56 * mrSges(4,1) + t55 * mrSges(4,2);
t87 = m(4) * pkin(2) + m(5) * t41 + mrSges(3,1) + t65 - t66;
t79 = pkin(5) + qJ(3);
t86 = -m(4) * qJ(3) - m(5) * t79 + mrSges(3,2) - mrSges(4,3) - mrSges(5,3);
t58 = sin(qJ(4));
t59 = cos(qJ(4));
t30 = t55 * t59 + t56 * t58;
t26 = t30 * qJD(2);
t84 = t26 / 0.2e1;
t83 = m(2) + m(3);
t82 = Ifges(5,4) * t26;
t78 = qJ(3) * qJD(2);
t32 = t55 * qJD(1) + t56 * t78;
t77 = qJDD(2) * t55;
t76 = qJDD(2) * t56;
t36 = t79 * t55;
t29 = -t55 * t58 + t56 * t59;
t27 = t29 * qJD(4);
t12 = qJD(2) * t27 + qJDD(2) * t30;
t28 = t30 * qJD(4);
t13 = -qJD(2) * t28 + qJDD(2) * t29;
t69 = -t13 * mrSges(5,1) + t12 * mrSges(5,2);
t54 = pkin(6) + qJ(2);
t48 = sin(t54);
t50 = cos(t54);
t68 = -g(1) * t48 + g(2) * t50;
t67 = -mrSges(4,1) * t76 + mrSges(4,2) * t77;
t46 = t56 * qJD(1);
t20 = -qJD(2) * t36 + t46;
t21 = pkin(5) * qJD(2) * t56 + t32;
t5 = t20 * t59 - t21 * t58;
t6 = t20 * t58 + t21 * t59;
t63 = -(-t55 * t78 + t46) * t55 + t32 * t56;
t37 = t79 * t56;
t14 = -t36 * t59 - t37 * t58;
t15 = -t36 * t58 + t37 * t59;
t42 = -qJDD(2) * pkin(2) + qJDD(3);
t35 = -qJD(2) * t41 + qJD(3);
t33 = -qJDD(2) * t41 + qJDD(3);
t25 = t29 * qJD(2);
t22 = Ifges(5,4) * t25;
t19 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t26;
t18 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t25;
t17 = pkin(5) * t76 + t24;
t16 = t44 + (-pkin(5) * qJDD(2) - t38) * t55;
t10 = Ifges(5,1) * t26 + Ifges(5,5) * qJD(4) + t22;
t9 = Ifges(5,2) * t25 + Ifges(5,6) * qJD(4) + t82;
t8 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t13;
t7 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t12;
t4 = -qJD(3) * t30 - qJD(4) * t15;
t3 = qJD(3) * t29 + qJD(4) * t14;
t2 = -qJD(4) * t6 + t16 * t59 - t17 * t58;
t1 = qJD(4) * t5 + t16 * t58 + t17 * t59;
t11 = [t27 * t18 - t28 * t19 + t29 * t7 + t30 * t8 + t83 * qJDD(1) + m(4) * (t23 * t56 + t24 * t55) + m(5) * (t1 * t30 + t2 * t29 + t27 * t6 - t28 * t5) + (-m(4) - m(5) - t83) * g(3); m(5) * (t1 * t15 + t14 * t2 + t3 * t6 - t33 * t41 + t4 * t5) + (mrSges(5,2) * t33 - mrSges(5,3) * t2 + Ifges(5,1) * t12 + Ifges(5,4) * t13 + Ifges(5,5) * qJDD(4)) * t30 + (-mrSges(5,1) * t33 + mrSges(5,3) * t1 + Ifges(5,4) * t12 + Ifges(5,2) * t13 + Ifges(5,6) * qJDD(4)) * t29 + m(4) * (-pkin(2) * t42 + t88 * qJ(3) + t63 * qJD(3)) + (t87 * t48 + t86 * t50) * g(1) + (t86 * t48 - t87 * t50) * g(2) + (-t27 * t5 - t28 * t6) * mrSges(5,3) + t35 * (mrSges(5,1) * t28 + mrSges(5,2) * t27) + t25 * (Ifges(5,4) * t27 - Ifges(5,2) * t28) / 0.2e1 + qJD(4) * (Ifges(5,5) * t27 - Ifges(5,6) * t28) / 0.2e1 + (Ifges(5,1) * t27 - Ifges(5,4) * t28) * t84 - t28 * t9 / 0.2e1 + t27 * t10 / 0.2e1 + t14 * t7 + t15 * t8 + t3 * t18 + t4 * t19 + Ifges(3,3) * qJDD(2) + (Ifges(4,4) * t55 + Ifges(4,2) * t56) * t76 + (Ifges(4,1) * t55 + Ifges(4,4) * t56) * t77 + t42 * t66 - pkin(2) * t67 - t41 * t69 + (t90 * t38 + t88) * mrSges(4,3); -t25 * t18 + t26 * t19 - t90 * qJD(2) ^ 2 * mrSges(4,3) + t67 + t69 + (-t6 * t25 + t5 * t26 + t33 + t68) * m(5) + (-qJD(2) * t63 + t42 + t68) * m(4); Ifges(5,5) * t12 + Ifges(5,6) * t13 + Ifges(5,3) * qJDD(4) - t1 * mrSges(5,2) + t2 * mrSges(5,1) - t35 * (mrSges(5,1) * t26 + mrSges(5,2) * t25) - t26 * (Ifges(5,1) * t25 - t82) / 0.2e1 + t9 * t84 - qJD(4) * (Ifges(5,5) * t25 - Ifges(5,6) * t26) / 0.2e1 - t5 * t18 + t6 * t19 - g(3) * t65 + (t25 * t5 + t26 * t6) * mrSges(5,3) - (-Ifges(5,2) * t26 + t10 + t22) * t25 / 0.2e1 + (g(1) * t50 + g(2) * t48) * (mrSges(5,1) * t47 + mrSges(5,2) * t49);];
tau = t11;
